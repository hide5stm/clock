#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>

const int N = 60; //粒子数
int n, q, p;

double a_large = 0.08,a_small = 0.01;//radius of disk
double m_large = 1.0,m_small = 1.0;
double h = 3.0;
double rate_large = 0.2;
//double a = 0.01;//粒子の半径
double e = 0.95,e_wall = 1.0;//それぞれ粒子同士、壁との反発係数
double g = 1.0;//規格化された重力加速度
double Xmin = -3.0,Xmax = 3.0,X0 = 0.0;//左右の壁の位置
double Ymin = 0.0,Ymax = 6.0;//底面とセルの最高点の位置
double U = 0.54;//床面の振動する速度
double V0 = 0.5;//初期条件での速度分布の標準偏差
int N_cell_x = 16,N_cell_y = 12;//x,y方向のセルの分割数
double T = 50.0;//シミュレーション終了時刻
double epsilon = 0.000001;//数値誤差の影響を除くために入れている
int step_gif = 400;//gifアニメーションのステップ数

//完全平衡木のノードの構造体
struct NODE {
	int number;//対応する粒子番号
	double time;//対応する粒子の予測された最短衝突時刻
	struct NODE *left;//上と左右をつなぐ
	struct NODE *right;
	struct NODE *parent;
};

//ある粒子のイベントの詳細(衝突時刻・相手)を記録
struct EVENT {
	double time;
	int number_particle;
};

//粒子に関する情報をまとめる
struct PARTICLE {
	double x,y,u,v;//位置・速度
	double tau;//粒子の固有時間を記録，Delayed State Algorithm(DSA)による高速化のために必要
	double next;//自分の次に同じセルに入った粒子の番号
	double a;
	double m;
	char size[10];
	struct EVENT event;//次に予定されるイベント
};

struct CELL {
	int first;//セルへの登録時に最も早く登録された粒子
	int last;
};

int intpow(int a, int b);
struct EVENT Predictions(struct PARTICLE pl[N], struct CELL cell[N_cell_x][N_cell_y], double t, int i);
void CBT_build(struct NODE *node[n+1][2 * p + 2 * q], struct PARTICLE pl[N]);
void CBT_update(struct NODE *entry[n+1][2 * p + 2 * q], double time_new, int i_new);
void status_initialize(struct PARTICLE pl[N]);
void cell_register(struct PARTICLE pl[N], struct CELL cell[N_cell_x][N_cell_y]);
int set(struct PARTICLE partile[N], int i);
double r_distance(struct PARTICLE particle1, struct PARTICLE particle2);
double v_distance(struct PARTICLE particle1, struct PARTICLE particle2);
double Uniform(void);
double rand_normal( double mu, double sigma );
int getcell_x(double x, double cell_length_x);//ellはcell_lengthのこと
int getcell_y(double y, double cell_length_y);
void Free_evolution(struct PARTICLE *particle, double t);
void G1(struct PARTICLE *particle, int j);
void G2(struct PARTICLE *particle1, struct PARTICLE *particle2);
double T_DWC(struct PARTICLE particle, double t, int j);
double T_DDC(struct PARTICLE particle1, struct PARTICLE particle2, double t);
double NextEvent(struct PARTICLE pl[N], struct CELL cell[N_cell_x][N_cell_y], struct NODE *node[n+1][2 * p + 2 * q], int i_current, int j_current);
double t_cell_update(struct PARTICLE particle, int j_current, double t_cell_old, double *v_max);
double Vmax(struct PARTICLE pl[N]);
void MaskUpdate(struct PARTICLE pl[N], struct CELL cell[N_cell_x][N_cell_y], struct NODE *node[n+1][2 * p + 2 * q], int i_current, double t);
double EEPGM(struct PARTICLE pl[N], struct CELL cell[N_cell_x][N_cell_y], struct NODE *node[n+1][2 * p + 2 * q], double t, double *v_max);
void SolveCubicEquation(double complex x[3], double a, double b, double c, double d);
double SolveQuarticEquation(double a, double b, double c, double d, double e);
double Cuberoot(double x);
double complex Squreroot(double complex x);
double intmax(int x, int y);
double m_cal(struct PARTICLE pl[N]);
double bias(struct PARTICLE pl[N],char c[10]);


int main(void) {
	n = (int)ceil(log2(N)); //完全平衡木の深さ
	q = (int)pow(2, n) - N;
	p = (N - q)/2;

	FILE *fp_large, *fp_small, *fp_bias, *fp_setting;
	if ((fp_large = fopen("large.txt", "w")) == NULL) {
		printf("fp_large open error\n");
	}
	if ((fp_small = fopen("small.txt", "w")) == NULL) {
		printf("fp_small open error\n");
	}
	if ((fp_bias = fopen("bias.txt", "w")) == NULL) {
		printf("fp_bias open error\n");
	}
	if ((fp_setting = fopen("setting.txt", "w"))==NULL) {//gif生成時に必要なパラメータを格納する
		printf("file open error\n");
	}

	int i_current, j_current;//現在注目している粒子のペア, j_current<0:壁, j_current>=0:粒子
	double v_max = 0.0;//最大速度を保存、セルの更新のために必要
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x, cell_length_y =(Ymax-Ymin)/(double)N_cell_y;//1つのセルの高さと横幅
	double t=0.0, dt=0.01, trec=0.0, dtrec = (double)T/(double)step_gif;//dtrecはgif生成の時間間隔
	double t_cell=0.0, t_cell_old=0.0;//セルの更新時刻
	double bias_large, bias_small;
	int N_large = 0, N_small = 0;
	//gif生成に必要な情報を出力
	fprintf(fp_setting, "Xmin = %lf\nXmax = %lf\nYmin = %lf\nYmax = %lf\nh = %lf\n", Xmin, Xmax, Ymin, Ymax, h);
	fprintf(fp_setting, "n1 = %d\ndt = %lf\n", step_gif, dtrec);
	fprintf(fp_setting, "file = \"position(U=%lf)\"\n", U);

	srand((unsigned) time(NULL));
	struct PARTICLE pl[N];//粒子を表す構造体
	struct CELL cell[N_cell_x][N_cell_y];
	status_initialize(pl);//位置や速度の初期化
	v_max = Vmax(pl);//粒子の最大速度
	t_cell = (cell_length_y-2.0*a_large)/(2.0*v_max);//この時間までにマスク外からの衝突はありえない
	cell_register(pl, cell);//粒子をセルに登録する, nextofの初期化

	for (int i = 0; i < N; i++) {
		if (strcmp(pl[i].size, "small") == 0) {
			//printf("small\n");
			N_small += 1;
		} else {
			N_large += 1;
			//printf("large\n");
		}
	}

	printf("n=%d p=%d q=%d\n", n, p, q);
	struct NODE *node[n+1][2 * p + 2 * q];//完全平衡木(あるいはトーナメント)を表す構造体
	//nodeの実体化、もうちょっとおしゃれにしたい
	for (int i = 0; i <= n; i++) {
		for (int j = 0; j < 2 * p + 2 * q; j++) {
			node[i][j] = (struct NODE *)malloc(sizeof(struct NODE));//領域を確保する
		}
	}

	for (int i = 0; i < N; i++) {
		pl[i].event = Predictions(pl, cell, t, i);//それぞれの粒子の次のイベントを予測
	}
	CBT_build(node, pl);//Complete Binary Treeを組み立てる
	printf("set up ok\n");
	printf("N_large:%d, N_small:%d\n", N_large, N_small);

	while(t <= T) {
		//NEXT EVENTの検索
		i_current = node[0][0]->number;//決勝のノードは最短の時間で衝突する粒子を示す
		j_current = pl[i_current].event.number_particle;//i_currentの衝突相手(これは壁の可能性もある)
		t = NextEvent(pl, cell, node, i_current, j_current);//NEXT EVENTを処理しtとparticle, cell, nodeを更新
		t_cell = t_cell_update(pl[i_current], j_current, t_cell_old, &v_max);//t_cellとv_maxの更新
		//i_current, j_currentと衝突する予定だった粒子がいた場合はその粒子のeventはinvalidになってしまうので新しくeventを作る
		//そのような粒子は同じマスク内にしか存在しないはずなのでその中で探索
		MaskUpdate(pl, cell, node, i_current, t);//i_currentの周りの粒子でinvalidなものがあればアップデート
		if (j_current >= 0) {//jについても同様
			MaskUpdate(pl, cell, node, j_current, t);
		}

		//EEPGM マスク外の粒子とも衝突する可能性が生じるので登録し直す
		if (t >= t_cell) {
			t_cell_old = t;
			t_cell = EEPGM(pl, cell, node, t, &v_max);
			//床に粒子がめり込んでいたらこのエラーが生じる
//			for (int i = 0; i < N; i++) {
//				if (pl[i].y < Ymin+a-epsilon) {
//					printf("i=%d:error\n", i);
//					printf("%lf %lf %lf %lf\n", pl[i].x, pl[i].y, pl[i].u, pl[i].v);
//					printf("%lf %d\n", pl[i].event.time, pl[i].event.number_particle);
//					G1(&pl[i], -3);
//					pl[i].event = Predictions(pl, cell, t, i);
//					CBT_update(node, pl[i].event.time, i);
//					MaskUpdate(pl, cell, node, i, t);
//				}
//			}

		}
		//粒子の位置の出力
		if ((t > trec)&&(t < T)) {
			t_cell_old = t;
			t_cell = EEPGM(pl, cell, node, t, &v_max);
			for (int i = 0; i < N; i++) {
				if (strcmp(pl[i].size, "small") == 0) {
					fprintf(fp_small, "%lf %lf\n", pl[i].x, pl[i].y);
				} else {
					fprintf(fp_large, "%lf %lf\n", pl[i].x, pl[i].y);
				}
			}
			fprintf(fp_small, "\n\n");
			fprintf(fp_large, "\n\n");
			bias_large = bias(pl, "large")/(double)N_large;
			bias_small = bias(pl, "small")/(double)N_small;
			fprintf(fp_bias, "%lf %lf %lf\n", t, bias_large, bias_small);
			printf("t = %lf, bias_large = %lf, bias_small = %lf\n", t, bias_large, bias_small);
			trec += dtrec;

		}

	}

	fclose(fp_large);
	fclose(fp_small);
	fclose(fp_bias);
	fclose(fp_setting);
	system("gnuplot gif.txt\n");
	system("gnuplot bias_plt.txt\n");
	return 0;
}

//pow()の整数版, トーナメントを作成するときに必要
int intpow(int a, int b) {
	return (int)pow(a, b);
}

//粒子iに対して最短で生じるイベントを予測して出力
struct EVENT Predictions(struct PARTICLE pl[N],
						 struct CELL cell[N_cell_x][N_cell_y],
						 double t, int i) {
	double t_min = 2.0*T, t_temp;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x, cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int j_col, j;
	struct PARTICLE pl_j;
	struct EVENT L;

	for (j = -5; j < 0; j++) {//壁との衝突時間を確認
		t_temp = T_DWC(pl[i], pl[i].tau, j);
		if ((t_temp > t) & ( t_temp < t_min)) {
			t_min = t_temp;
			j_col = j;
		}
	}

	int cell_x = getcell_x(pl[i].x, cell_length_x), cell_y = getcell_y(pl[i].y, cell_length_y);
	for (int c1 = -1; c1 <= 1; c1++) {//近辺の粒子との衝突時間を確認
		for (int c2 = -1; c2 <= 1; c2++) {
			if ((((cell_x+c1 >= 0) && (cell_x+c1 < N_cell_x)) && (cell_y+c2 >= 0)) && (cell_y+c2 < N_cell_y)) {
				j = cell[cell_x+c1][cell_y+c2].first;//リンクリスト構造を利用して粒子を捜索
				while(j >= 0) {
					pl_j = pl[j];
					Free_evolution(&pl_j, pl[i].tau-pl[j].tau);//相手粒子jの時間をiとそろえる
					t_temp = T_DDC(pl[i], pl_j, pl[i].tau);//衝突にかかる時間を計算
					//現在時刻より遅く，t_minよりも早い時間であればt_minの更新
					if ((t_temp > t) && (t_temp < t_min)) {
						t_min = t_temp;
						j_col = j;//そのときの相手jも記録
					}
					j = pl[j].next;
				}
			}
		}
	}
	L.time = t_min;
	L.number_particle = j_col;
	//時間と相手の情報を出力
	return L;
}

//CBTを作る
void CBT_build(struct NODE *node[n+1][2 * p + 2 * q],
			   struct PARTICLE pl[N]) {
	//initialization for bottom nodes
	for (int i = 0; i < 2 * p + 2 * q; i++) {
		if (i < 2*p) {
			node[n][i]->time = pl[i].event.time;
			node[n][i]->number = i;
		} else {
			node[n][i]->time = pl[2*p+(i-2*p)/2].event.time;
			node[n][i]->number = 2*p+(i-2*p)/2;
		}
	}
	//tournament
	for (int j = n - 1; j >= 0; j--) {
		for (int i = 0; i <= intpow(2, j) - 1; i++) {
			node[j][i]->left = node[j + 1][2 * i];
			node[j + 1][2 * i]->parent = node[j][i];
			node[j][i]->right = node[j + 1][2 * i + 1];
			node[j + 1][2 * i + 1]->parent = node[j][i];
			if (node[j + 1][2 * i]->time <= node[j + 1][2 * i + 1]->time) {
				node[j][i]->time = node[j + 1][2 * i]->time;
				node[j][i]->number = node[j + 1][2 * i]->number;
			} else {
				node[j][i]->time = node[j + 1][2 * i + 1]->time;
				node[j][i]->number = node[j + 1][2 * i + 1]->number;
			}
		}
	}
}

//i_newの情報を更新する
void CBT_update(struct NODE *entry[n+1][2 * p + 2 * q],
				double time_new, int i_new) {
	struct NODE *entry_now, hoge_now;
	if (i_new < 2*p) {
		entry[n][i_new]->time = time_new;
		entry_now = entry[n][i_new];
	} else {
		entry[n][2*i_new-2*p]->time = time_new;
		entry[n][2*i_new-2*p+1]->time = time_new;
		entry_now = entry[n][2*i_new-2*p];
	}
	while(entry_now->parent != NULL) {
		entry_now = entry_now->parent;
		if (entry_now->left->time < entry_now->right->time) {
			entry_now->time = entry_now->left->time;
			entry_now->number = entry_now->left->number;
		} else {
			entry_now->time = entry_now->right->time;
			entry_now->number = entry_now->right->number;
		}
	}
}

//粒子の初期条件を決める
void status_initialize(struct PARTICLE pl[N]) {
	double prob;
	int i;
	for (i = 0; i < N; i++) {
		prob = Uniform();
		if (prob < rate_large) {
			pl[i].m = m_large;
			pl[i].a = a_large;
			strcpy(pl[i].size, "large");
			prob = Uniform();
			pl[i].x = (X0+pl[i].a)*prob+(Xmax-pl[i].a)*(1-prob);
			prob = Uniform();
			pl[i].y = (Ymin+pl[i].a)*prob+(Ymax-pl[i].a)*(1-prob);
			while(set(pl, i) == 0) {//もし重なっている粒子があったときは重なりがなくなるまで登録し直す
				prob = Uniform();
				pl[i].x = (X0+pl[i].a)*prob+(Xmax-pl[i].a)*(1-prob);
				prob = Uniform();
				pl[i].y = (Ymin+pl[i].a)*prob+(Ymax-pl[i].a)*(1-prob);
			}
		} else {
			pl[i].m = m_small;
			pl[i].a = a_small;
			strcpy(pl[i].size, "small");
			prob = Uniform();
			pl[i].x = (Xmin+pl[i].a)*prob+(X0-pl[i].a)*(1-prob);
			prob = Uniform();
			pl[i].y = (Ymin+pl[i].a)*prob+(Ymax-pl[i].a)*(1-prob);
			while(set(pl, i) == 0) {//もし重なっている粒子があったときは重なりがなくなるまで登録し直す
				prob = Uniform();
				pl[i].x = (Xmin+pl[i].a)*prob+(X0-pl[i].a)*(1-prob);
				prob = Uniform();
				pl[i].y = (Ymin+pl[i].a)*prob+(Ymax-pl[i].a)*(1-prob);
			}
		}
		/*
		prob = Uniform();
		pl[i].x = (Xmin+pl[i].a)*prob+(Xmax-pl[i].a)*(1-prob);
		prob = Uniform();
		pl[i].y = (Ymin+pl[i].a)*prob+(0.5*Ymax-pl[i].a)*(1-prob);

		while(set(pl, i) == 0) {//もし重なっている粒子があったときは重なりがなくなるまで登録し直す
			prob = Uniform();
			pl[i].x = (Xmin+pl[i].a)*prob+(Xmax-pl[i].a)*(1-prob);
			prob = Uniform();
			pl[i].y = (Ymin+pl[i].a)*prob+(Ymax-pl[i].a)*(1-prob);
		}
		*/
		pl[i].u = rand_normal(0.0, V0);
		pl[i].v = rand_normal(0.0, V0);
		pl[i].next = -1;
		pl[i].tau = 0.0;
		pl[i].event.time = 2.0*T;
		pl[i].event.number_particle = -1;
	}
}

//すべての粒子をセルに登録し直す
void cell_register(struct PARTICLE pl[N],
				   struct CELL cell[N_cell_x][N_cell_y]) {
	int i, x, y, lastPrev;
	//initialize pl.next and cell
	for (i = 0; i < N; i++) {
		pl[i].next = -1;
	}
	for (x = 0; x < N_cell_x; x++) {
		for (y = 0; y < N_cell_y; y++) {
			cell[x][y].first = -1;
			cell[x][y].last = -1;
		}
	}
	//リンクリスト構造の作成
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x;
	double cell_length_y = (Ymax-Ymin)/(double)N_cell_y;

	for (i = 0; i < N; i++) {
		x = getcell_x(pl[i].x, cell_length_x);
		y = getcell_y(pl[i].y, cell_length_y);
		lastPrev = cell[x][y].last;
		cell[x][y].last = i;

		if (lastPrev == -1) {
			cell[x][y].first = i;
		} else {
			pl[lastPrev].next = i;
		}
	}
}

//重なりなく粒子を置くことに成功していれば1, 失敗していれば0を返す
int set(struct PARTICLE pl[N], int i) {
	int j, r = 1;
	double d;

	if (fabs(pl[i].x) < pl[i].a) {
		r = 0;
	}
	for (j = 0; j <= i - 1; j++) {
		d = r_distance(pl[i], pl[j]);
		if (d <= pl[i].a + pl[j].a) {
			r = 0;
			break;
		}
	}
	return r;
}

double r_distance(struct PARTICLE pl1, struct PARTICLE pl2) {//2つの粒子の距離を計算
	double d;
	d = sqrt(pow(pl1.x-pl2.x, 2.0)+pow(pl1.y-pl2.y, 2.0));
	return d;
}

double v_distance(struct PARTICLE pl1, struct PARTICLE pl2) {//2つの粒子の速度ベクトルの差の大きさを計算
	double d;
	d = sqrt(pow(pl1.u-pl2.u, 2.0)+pow(pl1.v-pl2.v, 2.0));
	return d;
}

double Uniform(void) {//0から1の一様乱数を生成
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_normal( double mu, double sigma ) {//平均mu, 標準偏差sigmaの正規分布を生成
	double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
}

int getcell_x(double x, double cell_length_x) {//xという位置のセルの番号を返す
	if ((x < Xmin+a_small)||(Xmax-a_small < x)) {
		printf("x is out of range\n");
	}
	return (int)((x-Xmin)/cell_length_x);
}

int getcell_y(double y, double cell_length_y) {//yという位置のセルの番号を返す
	if (y < Ymin) {
		printf("error:y<0(%lf)\n", y);
		return 0;
	} else if (y>Ymax) {
		return N_cell_y-1;//Ymaxよりも高い位置の粒子は一番高いセルに登録
	} else {
		return (int)(y/cell_length_y);
	}
}

//ある粒子を時間tだけ時間発展させる
void Free_evolution(struct PARTICLE *pl, double t) {
	pl->x += (pl->u)*t;
	pl->y += (pl->v)*t-0.5*g*t*t;
	pl->v += -g*t;
	pl->tau += t;//固有時間の更新も必要なことに注意
}

//粒子と壁の衝突処理を行う
void G1(struct PARTICLE *pl, int j) {
	double temp;
	if ((j == -1) || (j == -2)) {//collision with R or L wall
		pl->u = -e_wall*pl->u;
		if (j == -1) {
			pl->x = Xmax-pl->a-epsilon;//このepsilon処理はgetcell_xのときなどに必要になる
		} else {
			pl->x = Xmin+pl->a+epsilon;//ここから修正すること
		}
	} else if (j == -3) {//collision with Bottom wall
		pl->v = (1+e_wall)*U-e_wall*pl->v;
		pl->y = Ymin+pl->a+epsilon;
	} else if (j == -4) {//collision with Center Wall
		pl->x = (pl->x/fabs(pl->x))*pl->a;
		if (pl->y < h) {
		//if (fabs(pl->y-h) > ell) {//really collision
			pl->u = -e_wall*pl->u;
		}
	} else if (j == -5) {//collision with upper edge
		double nx = (pl->x-X0)/pl->a, ny = (pl->y-h)/pl->a;
		double utemp = pl->u, vtemp = pl->v;
		pl->u = -2.0*(utemp*nx+vtemp*ny)*nx+utemp;
		pl->v = -2.0*(utemp*nx+vtemp*ny)*ny+vtemp;
		if (nx*nx+ny*ny < 1.0-epsilon) {//ふちにめり込んでいるときにこのエラーが生じる
			printf("check:%lf\n", nx*nx+ny*ny);
			printf("upper:%lf %lf\n", pl->x, pl->y);
			printf("velocity:%lf %lf -> %lf %lf\n", utemp, vtemp, pl->u, pl->v);
			printf("energy:%lf -> %lf\n", utemp*utemp+vtemp*vtemp, pl->u*pl->u+pl->v*pl->v);
		}
	}
}

//粒子同士の衝突処理
void G2(struct PARTICLE *pl1, struct PARTICLE *pl2) {
	double d = r_distance(*pl1, *pl2);
	double Utemp1 = pl1->u;
	double Vtemp1 = pl1->v;
	double Utemp2 = pl2->u;
	double Vtemp2 = pl2->v;
	double Cx = (pl1->x-pl2->x)/d;
	double Cy = (pl1->y-pl2->y)/d;

	pl1->u = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cx+Utemp1;
	pl1->v = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cy+Vtemp1;
	pl2->u = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cx+Utemp2;
	pl2->v = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cy+Vtemp2;
}

//粒子同士の衝突(DDC)の時間を計算
double T_DDC(struct PARTICLE pl1, struct PARTICLE pl2, double t) {
	double r_relative, v_relative, b, hoge;
	double tau = t;
	double x1 = pl1.x, x2 = pl2.x, y1 = pl1.y, y2 = pl2.y;
	double u1 = pl1.u, u2 = pl2.u, v1 = pl1.v, v2 = pl2.v;
	r_relative = r_distance(pl1, pl2);
	v_relative = v_distance(pl1, pl2);
	b = (x1-x2)*(u1-u2)+(y1-y2)*(v1-v2);
	hoge = b*b-v_relative*v_relative*(r_relative*r_relative-pow(pl1.a+pl2.a, 2.0));
	if (hoge > 0.0) {
		tau += -(b+sqrt(hoge))/(v_relative*v_relative);
	} else {
		tau += T;
	}
	return tau;
}

//粒子と壁の衝突の時間を計算
double T_DWC(struct PARTICLE pl, double t, int j) {
	double tau = t;
	if (j==-1) {//collision with RIGHT wall(-1)
		if (pl.u > 0.0) {
			tau += (Xmax-pl.a-pl.x)/pl.u;
		} else {
			tau += 2.0*T;
		}
	} else if (j == -2) {//collision with LEFT wall(-2)
		if (pl.u<0.0) {
			tau += (Xmin+pl.a-pl.x)/pl.u;
		} else {
			tau += 2.0*T;
		}
	} else if (j == -3) {//collision with BOTTOM wall(-3)
		tau += (pl.v+sqrt(pl.v*pl.v+2*g*(pl.y-pl.a)))/g;
	} else if (j==-4) {
		if (pl.x > pl.a) {
			tau += -(pl.x-pl.a)/pl.u;
		} else if (pl.x < -pl.a) {
			tau += (-pl.x-pl.a)/pl.u;
		}
	} else {//j=-5
		if (fabs(pl.x-X0) <= pl.a) {
			tau += SolveQuarticEquation(g * g /4.0,
										-pl.v * g,
										-(pl.y - h) * g + pl.u * pl.u + pl.v * pl.v,
										2.0 * ((pl.x - X0) * pl.u + (pl.y - h) * pl.v),
										pow(pl.x-X0, 2.0) + pow(pl.y-h, 2.0) - pl.a * pl.a);
		}
	}
	if (tau < t) {
		return 2.0 * T;
	} else {
		return tau;
	}
}

//i_currentとj_currentのイベントを実際に行い，その時刻を返す関数
double NextEvent(struct PARTICLE pl[N],
				 struct CELL cell[N_cell_x][N_cell_y],
				 struct NODE *node[n+1][2 * p + 2 * q],
				 int i_current, int j_current) {
	double t = pl[i_current].event.time;
	Free_evolution(&pl[i_current], t-pl[i_current].tau);//i_currentの時間発展
	if (j_current >= 0) {//Disk Disk Collision
		Free_evolution(&pl[j_current], t - pl[j_current].tau);//j_currentの時間発展
		G2(&pl[i_current], &pl[j_current]);//粒子同士の衝突処理
		if (r_distance(pl[i_current], pl[j_current]) < pl[i_current].a + pl[j_current].a - epsilon) {
			printf("%d %d is too close!!\n", i_current, j_current);
			printf("distance = %lf\n", r_distance(pl[i_current], pl[j_current]));
		}
	}
	if (j_current < 0) {//Disk Wall Collision
		G1(&pl[i_current], j_current);//壁との衝突処理
	}
	pl[i_current].event = Predictions(pl, cell, t, i_current);//i_currentのイベント更新
	CBT_update(node, pl[i_current].event.time, i_current);//i_currentのnodeアップデート
	if (j_current >= 0) {//j_currentについても同様
		pl[j_current].event = Predictions(pl, cell, t, j_current);
		CBT_update(node, pl[j_current].event.time, j_current);
	}
	return t;
}

//セルの更新時刻を計算
double t_cell_update(struct PARTICLE pl,
					 int j_current, double t_cell_old, double *v_max) {
	double t_cell, dt_cell;
	double cell_length_x = (Xmax-Xmin) / (double)N_cell_x;
	double cell_length_y = (Ymax-Ymin) / (double)N_cell_y;

	if (j_current == -3) {
		if (*v_max * (*v_max) < pow(pl.u, 2.0) + pow(pl.v, 2.0)) {
			*v_max = sqrt(pow(pl.u, 2.0) + pow(pl.v, 2.0));
		}
	}
	dt_cell = (cell_length_y-2.0 * a_large)/(2.0 * (*v_max));
	t_cell = t_cell_old + dt_cell;
	return t_cell;
}

//最大速度を計算
double Vmax(struct PARTICLE pl[N]) {
	double v_max = 0.0;
	for (int i = 0; i < N; i++) {
		if (v_max*v_max < pl[i].u*pl[i].u+pl[i].v*pl[i].v) {
			v_max = sqrt(pl[i].u*pl[i].u+pl[i].v*pl[i].v);
		}
	}
	return v_max;
}

//同じマスクに含まれる粒子のイベントを更新
void MaskUpdate(struct PARTICLE pl[N],
				struct CELL cell[N_cell_x][N_cell_y],
				struct NODE *node[n+1][2 * p + 2 * q],
				int i_current, double t) {

	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x, cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int cell_x = getcell_x(pl[i_current].x, cell_length_x) , cell_y = getcell_y(pl[i_current].y, cell_length_y), j;
	for (int c1=-1;c1<=1;c1++) {
		for (int c2=-1;c2<=1;c2++) {
			if ((((cell_x+c1 >= 0) && (cell_x+c1 < N_cell_x)) && (cell_y+c2 >= 0)) && (cell_y+c2 < N_cell_y)) {
				j = cell[cell_x+c1][cell_y+c2].first;
				while(j >= 0) {
					if (pl[j].event.number_particle == i_current) {
						pl[j].event = Predictions(pl, cell, t, j);
						CBT_update(node, pl[j].event.time, j);
					}
					j = pl[j].next;
				}
			}
		}
	}
}

//すべての粒子に関して時間発展させたのちセルに登録し直す
double EEPGM(struct PARTICLE pl[N],
			 struct CELL cell[N_cell_x][N_cell_y],
			 struct NODE *node[n+1][2 * p + 2 * q],
			 double t, double *v_max) {

	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x, cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	double dt_cell, t_cell;
	double t_free, t_register, t_prediction, t_build, t_vmax;
	clock_t  start, end;
	start = clock();
	for (int i = 0; i < N; i++) {//現在時刻まで時間発展
		Free_evolution(&pl[i], t-pl[i].tau);
	}
	end = clock();
	t_free = (double)(end-start);

	start = clock();
	cell_register(pl, cell);//全粒子をセルに登録し直す
	end = clock();
	t_register = (double)(end-start);

	start = clock();
	for (int i = 0; i < N; i++) {//全粒子についてeventを計算し直す
		pl[i].event = Predictions(pl, cell, t, i);
	}
	end = clock();
	t_prediction = (double)(end-start);

	start = clock();
	CBT_build(node, pl);//CBTも最初から構成
	end = clock();
	t_build = (double)(end-start);

	start = clock();
	*v_max = Vmax(pl);
	end = clock();
	t_vmax = (double)(end-start);
	printf("free:%lf, register:%lf, prediction:%lf, build:%lf, vmax:%lf\n", t_free, t_register, t_prediction, t_build, t_vmax);
	dt_cell = (cell_length_y-2.0*a_large)/(2.0*(*v_max));
	t_cell = t+dt_cell;
	return t_cell;
}

void SolveCubicEquation(double complex x[3], double a, double b, double c, double d) {
	if (a == 0.0) {
		printf("Error:a = 0.0\n");
		printf("This equation is NOT Cubic.\n");
	}
	double A, B, C, p_temp, q_temp, D;
	A = b/a;
	B = c/a;
	C = d/a;
	p_temp = B-A*A/3.0;
	q_temp = 2.0*A*A*A/27.0-A*B/3.0+C;
	D = q_temp*q_temp/4.0+p_temp*p_temp*p_temp/27.0;
	if (D < 0.0) {//three real solutions
		double theta = atan2(sqrt(-D), -q_temp*0.5);
		x[0] = 2.0*sqrt(-p_temp/3.0)*cos(theta/3.0)-A/3.0;
		x[1] = 2.0*sqrt(-p_temp/3.0)*cos((theta+2.0*M_PI)/3.0)-A/3.0;
		x[2] = 2.0*sqrt(-p_temp/3.0)*cos((theta+4.0*M_PI)/3.0)-A/3.0;
	} else {//single real solution and two imaginary solutions(c.c)
		double u = Cuberoot(-q_temp*0.5+sqrt(D)), v = Cuberoot(-q_temp*0.5-sqrt(D));
		x[0] = u+v-A/3.0;
		x[1] = -0.5*(u+v)+sqrt(3.0)*0.5*1i*(u-v)-A/3.0;
		x[2] = -0.5*(u+v)-sqrt(3.0)*0.5*1i*(u-v)-A/3.0;
	}
}

double SolveQuarticEquation(double a, double b, double c, double d, double e) {
	if (a == 0.0) {
		printf("Error:a = 0.0\n");
		printf("This equation is NOT Quartic.\n");
	}

	double A, B, C, D, p_temp, q_temp, r_temp;
	A = b/a;
	B = c/a;
	C = d/a;
	D = e/a;
	p_temp = -6.0*pow(A/4.0, 2.0)+B;
	q_temp = 8.0*pow(A/4.0, 3.0)-2.0*B*A/4.0+C;
	r_temp = -3.0*pow(A/4.0, 4.0)+B*pow(A/4.0, 2.0)-C*A/4.0+D;
	double complex t_temp[3];
	SolveCubicEquation(t_temp, 1.0, -p_temp, -4.0*r_temp, 4.0*p_temp*r_temp-q_temp*q_temp);
	double t = creal(t_temp[0]);
	double complex m_temp = Squreroot(t-p_temp);//mは誤差なく実数or純虚数
	double tau = 2.0*T;
	double complex x[4];
	x[0] = (-m_temp+Squreroot(-t-p_temp+2.0*q_temp/m_temp))*0.5-A/4.0;
	x[1] = (-m_temp-Squreroot(-t-p_temp+2.0*q_temp/m_temp))*0.5-A/4.0;
	x[2] = (m_temp+Squreroot(-t-p_temp-2.0*q_temp/m_temp))*0.5-A/4.0;
	x[3] = (m_temp-Squreroot(-t-p_temp-2.0*q_temp/m_temp))*0.5-A/4.0;
	for (int i=0;i<4;i++) {
		if (cimag(x[i]) == 0.0) {//解が実
			if ((creal(x[i]) < tau) && (creal(x[i]) > 0.0)) {
				tau = creal(x[i]);
			}
		}
	}
	if (fabs(a*pow(tau, 4.0)+b*pow(tau, 3.0)+c*pow(tau, 2.0)+d*tau+e) > epsilon) {
		tau = 2.0*T;
	}
	return tau;
}

double Cuberoot(double x) {
	if (x > 0.0) {
		return pow(x, 1.0/3.0);
	} else {
		return -pow(-x, 1.0/3.0);
	}
}

double complex Squreroot(double complex x) {
	double complex y;
	double r = sqrt(creal(x)*creal(x)+cimag(x)*cimag(x));
	double theta = atan2(cimag(x), creal(x));
	//t-pは実数だからかならずこちらで処理される
	if (cimag(x) == 0.0) {
		if (creal(x) > 0.0) {
			y = sqrt(r);
		} else {
			y = sqrt(r)*1i;
		}
	} else {
		if (theta < 0.0) {
			theta += 2.0*M_PI;
		}
		double complex y = sqrt(r)*(cos(theta*0.5)+1i*sin(theta*0.5));
	}
	return y;
}


double m_cal(struct PARTICLE pl[N]) {
	int N_left = 0, N_right;
	for (int i = 0; i < N; i++) {
		if (pl[i].x < X0) {
			N_left += 1;
		}
	}
	N_right = N-N_left;
	double m = (intmax(N_left, N_right)-0.5*N)/(double)N;
	return m;
}

double intmax(int x, int y) {
	if (x > y) {
		return (double)x;
	} else {
		return (double)y;
	}
}

double bias(struct PARTICLE pl[N], char c[10]) {
	double temp=0.0;
	if (strcmp(c, "small")==0) {//small
		for (int i = 0; i < N; i++) {
			if (( strcmp(pl[i].size, "small")==0) && (pl[i].x < X0)) {
				temp += 1.0;
			}
		}
	} else {//large
		for (int i = 0; i < N; i++) {
			if (( strcmp(pl[i].size, "large")==0) && (pl[i].x < X0)) {
				temp += 1.0;
			}
		}
	}
	return temp;
}
