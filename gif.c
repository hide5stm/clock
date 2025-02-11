#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <complex.h>

#define N 600//粒子数
#define n ((int)ceil(log2(N)))//完全平衡木の深さ
#define q ((int)pow(2,n)-N)
#define p ((N-q)/2)
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
double T = 2000.0;//シミュレーション終了時刻
double epsilon = 0.000001;//数値誤差の影響を除くために入れている
int step_gif = 400;//gifアニメーションのステップ数

struct NODE{//完全平衡木のノードの構造体
	int number;//対応する粒子番号
	double time;//対応する粒子の予測された最短衝突時刻
	struct NODE *left;//上と左右をつなぐ
	struct NODE *right;
	struct NODE *parent;
};

struct EVENT{//ある粒子のイベントの詳細(衝突時刻・相手)を記録
	double time;
	int number_particle;
};

struct PARTICLE{//粒子に関する情報をまとめる
	double x,y,u,v;//位置・速度
	double tau;//粒子の固有時間を記録，Delayed State Algorithm(DSA)による高速化のために必要
	double next;//自分の次に同じセルに入った粒子の番号
	double a;
	double m;
	char size[10];
	struct EVENT event;//次に予定されるイベント
};

struct CELL{
	int first;//セルへの登録時に最も早く登録された粒子
	int last;
};

int intpow(int a,int b);
struct EVENT Predictions(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],double t,int i);
void CBT_build(struct NODE *node[n+1][2*p+2*q],struct PARTICLE particle[N]);
void CBT_update(struct NODE *entry[n+1][2*p+2*q],double time_new,int i_new);
void status_initialize(struct PARTICLE particle[N]);
void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]);
int set(struct PARTICLE partile[N],int i);
double r_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double v_distance(struct PARTICLE particle1,struct PARTICLE particle2);
double Uniform(void);
double rand_normal( double mu, double sigma );
int getcell_x(double x,double cell_length_x);//ellはcell_lengthのこと
int getcell_y(double y,double cell_length_y);
void Free_evolution(struct PARTICLE *particle,double t);
void G1(struct PARTICLE *particle,int j);
void G2(struct PARTICLE *particle1,struct PARTICLE *particle2);
double T_DWC(struct PARTICLE particle,double t,int j);
double T_DDC(struct PARTICLE particle1,struct PARTICLE particle2,double t);
double NextEvent(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,int j_current);
double t_cell_update(struct PARTICLE particle,int j_current,double t_cell_old,double *v_max);
double Vmax(struct PARTICLE particle[N]);
void MaskUpdate(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,double t);
double EEPGM(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],double t,double *v_max);
void SolveCubicEquation(double complex x[3],double a,double b,double c,double d);
double SolveQuarticEquation(double a,double b,double c,double d,double e);
double Cuberoot(double x);
double complex Squreroot(double complex x);
double intmax(int x,int y);
double m_cal(struct PARTICLE particle[N]);
double bias(struct PARTICLE particle[N],char c[10]);


int main(void){
	FILE *fp_large,*fp_small,*fp_bias,*fp_setting;
	if((fp_large = fopen("large.txt","w")) == NULL){
		printf("fp_large open error\n");
	}
	if((fp_small = fopen("small.txt","w")) == NULL){
		printf("fp_small open error\n");
	}
	if((fp_bias = fopen("bias.txt","w")) == NULL){
		printf("fp_bias open error\n");
	}
	if((fp_setting = fopen("setting.txt","w"))==NULL){//gif生成時に必要なパラメータを格納する
		printf("file open error\n");
	}
	
	int i_current,j_current;//現在注目している粒子のペア,j_current<0:壁,j_current>=0:粒子
	double v_max = 0.0;//最大速度を保存、セルの更新のために必要
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y =(Ymax-Ymin)/(double)N_cell_y;//1つのセルの高さと横幅
	double t=0.0,dt=0.01,trec=0.0,dtrec = (double)T/(double)step_gif;//dtrecはgif生成の時間間隔
	double t_cell=0.0,t_cell_old=0.0;//セルの更新時刻
	double bias_large,bias_small;
	int N_large = 0,N_small = 0;
	//gif生成に必要な情報を出力
	fprintf(fp_setting,"Xmin = %lf\nXmax = %lf\nYmin = %lf\nYmax = %lf\nh = %lf\n",Xmin,Xmax,Ymin,Ymax,h);
	fprintf(fp_setting,"n1 = %d\ndt = %lf\n",step_gif,dtrec);
	fprintf(fp_setting,"file = \"position(U=%lf)\"\n",U);
	
	srand((unsigned) time(NULL));
	struct PARTICLE particle[N];//粒子を表す構造体
	struct CELL cell[N_cell_x][N_cell_y];
	status_initialize(particle);//位置や速度の初期化
	v_max = Vmax(particle);//粒子の最大速度
	t_cell = (cell_length_y-2.0*a_large)/(2.0*v_max);//この時間までにマスク外からの衝突はありえない
	cell_register(particle,cell);//粒子をセルに登録する,nextofの初期化
	
	for(int i=0;i<N;i++){
		if(strcmp(particle[i].size,"small") == 0){
			//printf("small\n");
			N_small += 1;
		}else{
			N_large += 1;
			//printf("large\n");
		}
	}
	
	struct NODE *node[n+1][2*p+2*q];//完全平衡木(あるいはトーナメント)を表す構造体
	//nodeの実体化、もうちょっとおしゃれにしたい
	for(int i=0;i<=n;i++){
		for(int j=0;j<2*p+2*q;j++){
			node[i][j] = (struct NODE *)malloc(sizeof(struct NODE));//領域を確保する
		}
	}
	
	for(int i=0;i<N;i++){
		particle[i].event = Predictions(particle,cell,t,i);//それぞれの粒子の次のイベントを予測
	}
	CBT_build(node,particle);//Complete Binary Treeを組み立てる
	printf("set up ok\n");
	printf("N_large:%d,N_small:%d\n",N_large,N_small);
	
	while(t <= T){
		//NEXT EVENTの検索
		i_current = node[0][0]->number;//決勝のノードは最短の時間で衝突する粒子を示す
		j_current = particle[i_current].event.number_particle;//i_currentの衝突相手(これは壁の可能性もある)
		t = NextEvent(particle,cell,node,i_current,j_current);//NEXT EVENTを処理しtとparticle,cell,nodeを更新
		t_cell = t_cell_update(particle[i_current],j_current,t_cell_old,&v_max);//t_cellとv_maxの更新
		//i_current,j_currentと衝突する予定だった粒子がいた場合はその粒子のeventはinvalidになってしまうので新しくeventを作る
		//そのような粒子は同じマスク内にしか存在しないはずなのでその中で探索
		MaskUpdate(particle,cell,node,i_current,t);//i_currentの周りの粒子でinvalidなものがあればアップデート
		if(j_current >= 0){//jについても同様
			MaskUpdate(particle,cell,node,j_current,t);
		}
		
		//EEPGM マスク外の粒子とも衝突する可能性が生じるので登録し直す
		if(t >= t_cell){
			t_cell_old = t;
			t_cell = EEPGM(particle,cell,node,t,&v_max);
			//床に粒子がめり込んでいたらこのエラーが生じる
//			for(int i=0;i<N;i++){
//				if(particle[i].y < Ymin+a-epsilon){
//					printf("i=%d:error\n",i);
//					printf("%lf %lf %lf %lf\n",particle[i].x,particle[i].y,particle[i].u,particle[i].v);
//					printf("%lf %d\n",particle[i].event.time,particle[i].event.number_particle);
//					G1(&particle[i],-3);
//					particle[i].event = Predictions(particle,cell,t,i);
//					CBT_update(node,particle[i].event.time,i);
//					MaskUpdate(particle,cell,node,i,t);
//				}
//			}
			
		}
		//粒子の位置の出力
		if((t > trec)&&(t < T)){
			t_cell_old = t;
			t_cell = EEPGM(particle,cell,node,t,&v_max);
			for(int i=0;i<N;i++){
				if(strcmp(particle[i].size,"small") == 0){
					fprintf(fp_small,"%lf %lf\n",particle[i].x,particle[i].y);
				}else{
					fprintf(fp_large,"%lf %lf\n",particle[i].x,particle[i].y);
				}
			}
			fprintf(fp_small,"\n\n");
			fprintf(fp_large,"\n\n");
			bias_large = bias(particle,"large")/(double)N_large;
			bias_small = bias(particle,"small")/(double)N_small;
			fprintf(fp_bias,"%lf %lf %lf\n",t,bias_large,bias_small);
			printf("t = %lf, bias_large = %lf, bias_small = %lf\n",t,bias_large,bias_small);
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

int intpow(int a,int b){//pow()の整数版,トーナメントを作成するときに必要
	return (int)pow(a,b);
}

struct EVENT Predictions(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],double t,int i){//粒子iに対して最短で生じるイベントを予測して出力
	double t_min = 2.0*T,t_temp;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int j_col,j;
	struct PARTICLE particle_j;
	struct EVENT L;
	
	for(j=-5;j<0;j++){//壁との衝突時間を確認
		t_temp = T_DWC(particle[i],particle[i].tau,j);
		if((t_temp > t) && ( t_temp < t_min)){
			t_min = t_temp;
			j_col = j;
		}
	}
	
	int cell_x = getcell_x(particle[i].x,cell_length_x),cell_y = getcell_y(particle[i].y,cell_length_y);
	for(int c1=-1;c1<=1;c1++){//近辺の粒子との衝突時間を確認
		for(int c2=-1;c2<=1;c2++){
			if((((cell_x+c1 >= 0) && (cell_x+c1 < N_cell_x)) && (cell_y+c2 >= 0)) && (cell_y+c2 < N_cell_y)){
				j = cell[cell_x+c1][cell_y+c2].first;//リンクリスト構造を利用して粒子を捜索
				while(j >= 0){
					particle_j = particle[j];
					Free_evolution(&particle_j,particle[i].tau-particle[j].tau);//相手粒子jの時間をiとそろえる
					t_temp = T_DDC(particle[i],particle_j,particle[i].tau);//衝突にかかる時間を計算
					if((t_temp > t) && (t_temp < t_min)){//現在時刻より遅く，t_minよりも早い時間であればt_minの更新
						t_min = t_temp;
						j_col = j;//そのときの相手jも記録
					}
					j = particle[j].next;
				}
			}
		}
	}
	L.time = t_min;
	L.number_particle = j_col;
	return L;//時間と相手の情報を出力
}
void CBT_build(struct NODE *node[n+1][2*p+2*q],struct PARTICLE particle[N]){//CBTを作る
	int i,n_index;
	//initialization for bottom nodes
	for(i=0;i<2*p+2*q;i++){
		if(i < 2*p){
			node[n][i]->time = particle[i].event.time;
			node[n][i]->number = i;
		}else{
			node[n][i]->time = particle[2*p+(i-2*p)/2].event.time;
			node[n][i]->number = 2*p+(i-2*p)/2;
		}
	}
	//tournament
	for(n_index=n-1;n_index>=0;n_index--){
		for(i=0;i<=intpow(2,n_index)-1;i++){
			node[n_index][i]->left = node[n_index+1][2*i];
			node[n_index+1][2*i]->parent = node[n_index][i];
			node[n_index][i]->right = node[n_index+1][2*i+1];
			node[n_index+1][2*i+1]->parent = node[n_index][i];
			if(node[n_index+1][2*i]->time <= node[n_index+1][2*i+1]->time){
				node[n_index][i]->time = node[n_index+1][2*i]->time;
				node[n_index][i]->number = node[n_index+1][2*i]->number;
			}else{
				node[n_index][i]->time = node[n_index+1][2*i+1]->time;
				node[n_index][i]->number = node[n_index+1][2*i+1]->number;
			}
		}
	}
}

void CBT_update(struct NODE *entry[n+1][2*p+2*q],double time_new,int i_new){//i_newの情報を更新する
	struct NODE *entry_now,hoge_now;
	if(i_new < 2*p){
		entry[n][i_new]->time = time_new;
		entry_now = entry[n][i_new];
	}else{
		entry[n][2*i_new-2*p]->time = time_new;
		entry[n][2*i_new-2*p+1]->time = time_new;
		entry_now = entry[n][2*i_new-2*p];
	}
	while(entry_now->parent != NULL){
		entry_now = entry_now->parent;
		if(entry_now->left->time < entry_now->right->time){
			entry_now->time = entry_now->left->time;
			entry_now->number = entry_now->left->number;
		}else{
			entry_now->time = entry_now->right->time;
			entry_now->number = entry_now->right->number;
		}
	}
}

void status_initialize(struct PARTICLE particle[N]){//粒子の初期条件を決める
	double prob;
	int i;
	for(i=0;i<N;i++){
		prob = Uniform();
		if(prob < rate_large){
			particle[i].m = m_large;
			particle[i].a = a_large;
			strcpy(particle[i].size,"large");
			prob = Uniform();
			particle[i].x = (X0+particle[i].a)*prob+(Xmax-particle[i].a)*(1-prob);
			prob = Uniform();
			particle[i].y = (Ymin+particle[i].a)*prob+(Ymax-particle[i].a)*(1-prob);
			while(set(particle,i) == 0){//もし重なっている粒子があったときは重なりがなくなるまで登録し直す
				prob = Uniform();
				particle[i].x = (X0+particle[i].a)*prob+(Xmax-particle[i].a)*(1-prob);
				prob = Uniform();
				particle[i].y = (Ymin+particle[i].a)*prob+(Ymax-particle[i].a)*(1-prob);
			}
		}else{
			particle[i].m = m_small;
			particle[i].a = a_small;
			strcpy(particle[i].size,"small");
			prob = Uniform();
			particle[i].x = (Xmin+particle[i].a)*prob+(X0-particle[i].a)*(1-prob);
			prob = Uniform();
			particle[i].y = (Ymin+particle[i].a)*prob+(Ymax-particle[i].a)*(1-prob);
			while(set(particle,i) == 0){//もし重なっている粒子があったときは重なりがなくなるまで登録し直す
				prob = Uniform();
				particle[i].x = (Xmin+particle[i].a)*prob+(X0-particle[i].a)*(1-prob);
				prob = Uniform();
				particle[i].y = (Ymin+particle[i].a)*prob+(Ymax-particle[i].a)*(1-prob);
			}
		}
		/*
		prob = Uniform();
		particle[i].x = (Xmin+particle[i].a)*prob+(Xmax-particle[i].a)*(1-prob);
		prob = Uniform();
		particle[i].y = (Ymin+particle[i].a)*prob+(0.5*Ymax-particle[i].a)*(1-prob);
		
		while(set(particle,i) == 0){//もし重なっている粒子があったときは重なりがなくなるまで登録し直す
			prob = Uniform();
			particle[i].x = (Xmin+particle[i].a)*prob+(Xmax-particle[i].a)*(1-prob);
			prob = Uniform();
			particle[i].y = (Ymin+particle[i].a)*prob+(Ymax-particle[i].a)*(1-prob);
		}
		*/
		particle[i].u = rand_normal(0.0,V0);
		particle[i].v = rand_normal(0.0,V0);
		particle[i].next = -1;
		particle[i].tau = 0.0;
		particle[i].event.time = 2.0*T;
		particle[i].event.number_particle = -1;
	}
}

void cell_register(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y]){//すべての粒子をセルに登録し直す
	int i,cell_x,cell_y,lastPrev;
	//initialize particle.next and cell
	for(i=0;i<N;i++){
		particle[i].next = -1;
	}
	for(cell_x = 0;cell_x < N_cell_x;cell_x++){
		for(cell_y = 0;cell_y < N_cell_y;cell_y++){
			cell[cell_x][cell_y].first = -1;
			cell[cell_x][cell_y].last = -1;
		}
	}
	//リンクリスト構造の作成
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	for(i=0;i<N;i++){
		cell_x = getcell_x(particle[i].x,cell_length_x);
		cell_y = getcell_y(particle[i].y,cell_length_y);
		lastPrev = cell[cell_x][cell_y].last;
		cell[cell_x][cell_y].last = i;
		
		if(lastPrev == -1){
			cell[cell_x][cell_y].first = i;
		}else{
			particle[lastPrev].next = i;
		}
	}
}


int set(struct PARTICLE particle[N],int i){//重なりなく粒子を置くことに成功していれば1,失敗していれば0を返す
	int j,r=1;
	double d;
	
	if(fabs(particle[i].x) < particle[i].a){
		r = 0;
	}
	for(j=0;j<=i-1;j++){
		d = r_distance(particle[i],particle[j]);
		if(d <= particle[i].a+particle[j].a){
			r = 0;
			break;
		}
	}
	return r;
}

double r_distance(struct PARTICLE particle1,struct PARTICLE particle2){//2つの粒子の距離を計算
	double d;
	d = sqrt(pow(particle1.x-particle2.x,2.0)+pow(particle1.y-particle2.y,2.0));
	return d;
}

double v_distance(struct PARTICLE particle1,struct PARTICLE particle2){//2つの粒子の速度ベクトルの差の大きさを計算
	double d;
	d = sqrt(pow(particle1.u-particle2.u,2.0)+pow(particle1.v-particle2.v,2.0));
	return d;
}

double Uniform(void){//0から1の一様乱数を生成
	return ((double)rand()+1.0)/((double)RAND_MAX+2.0);
}

double rand_normal( double mu, double sigma ){//平均mu,標準偏差sigmaの正規分布を生成
	double z=sqrt( -2.0*log(Uniform()) ) * sin( 2.0*M_PI*Uniform() );
	return mu + sigma*z;
}

int getcell_x(double x,double cell_length_x){//xという位置のセルの番号を返す
	if((x < Xmin+a_small)||(Xmax-a_small < x)){
		printf("x is out of range\n");
	}
	return (int)((x-Xmin)/cell_length_x);
}

int getcell_y(double y,double cell_length_y){//yという位置のセルの番号を返す
	if(y < Ymin){
		printf("error:y<0(%lf)\n",y);
		return 0;
	}else if(y>Ymax){
		return N_cell_y-1;//Ymaxよりも高い位置の粒子は一番高いセルに登録
	}else{
		return (int)(y/cell_length_y);
	}
}

void Free_evolution(struct PARTICLE *particle,double t){//ある粒子を時間tだけ時間発展させる
	particle->x += (particle->u)*t;
	particle->y += (particle->v)*t-0.5*g*t*t;
	particle->v += -g*t;
	particle->tau += t;//固有時間の更新も必要なことに注意
}

void G1(struct PARTICLE *particle,int j){//粒子と壁の衝突処理を行う
	double temp;
	if((j == -1) || (j == -2)){//collision with R or L wall
		particle->u = -e_wall*particle->u;
		if(j == -1){
			particle->x = Xmax-particle->a-epsilon;//このepsilon処理はgetcell_xのときなどに必要になる
		}else{
			particle->x = Xmin+particle->a+epsilon;//ここから修正すること
		}
	}else if(j == -3){//collision with Bottom wall
		particle->v = (1+e_wall)*U-e_wall*particle->v;
		particle->y = Ymin+particle->a+epsilon;
	}else if(j == -4){//collision with Center Wall
		particle->x = (particle->x/fabs(particle->x))*particle->a;
		if(particle->y < h){
		//if(fabs(particle->y-h) > ell){//really collision
			particle->u = -e_wall*particle->u;
		}
	}else if(j == -5){//collision with upper edge
		double nx = (particle->x-X0)/particle->a,ny = (particle->y-h)/particle->a;
		double utemp = particle->u,vtemp = particle->v;
		particle->u = -2.0*(utemp*nx+vtemp*ny)*nx+utemp;
		particle->v = -2.0*(utemp*nx+vtemp*ny)*ny+vtemp;
		if(nx*nx+ny*ny < 1.0-epsilon){//ふちにめり込んでいるときにこのエラーが生じる
			printf("check:%lf\n",nx*nx+ny*ny);
			printf("upper:%lf %lf\n",particle->x,particle->y);
			printf("velocity:%lf %lf -> %lf %lf\n",utemp,vtemp,particle->u,particle->v);
			printf("energy:%lf -> %lf\n",utemp*utemp+vtemp*vtemp,particle->u*particle->u+particle->v*particle->v);
		}
	}
}


void G2(struct PARTICLE *particle1,struct PARTICLE *particle2){//粒子同士の衝突処理
	double d,Xtemp,Ytemp,Utemp1,Utemp2,Vtemp1,Vtemp2,Cx,Cy;
	d = r_distance(*particle1,*particle2);
	Utemp1 = particle1->u;
	Vtemp1 = particle1->v;
	Utemp2 = particle2->u;
	Vtemp2 = particle2->v;
	Cx = (particle1->x-particle2->x)/d;
	Cy = (particle1->y-particle2->y)/d;
	particle1->u = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cx+Utemp1;
	particle1->v = 0.5*(1+e)*((Utemp2-Utemp1)*Cx+(Vtemp2-Vtemp1)*Cy)*Cy+Vtemp1;
	particle2->u = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cx+Utemp2;
	particle2->v = 0.5*(1+e)*((Utemp1-Utemp2)*Cx+(Vtemp1-Vtemp2)*Cy)*Cy+Vtemp2;
}

double T_DDC(struct PARTICLE particle1,struct PARTICLE particle2,double t){//粒子同士の衝突(DDC)の時間を計算
	double r_relative,v_relative,b,hoge;
	double tau = t;
	double x1 = particle1.x ,x2 = particle2.x ,y1 = particle1.y , y2 = particle2.y;
	double u1 = particle1.u ,u2 = particle2.u ,v1 = particle1.v , v2 = particle2.v;
	r_relative = r_distance(particle1,particle2);
	v_relative = v_distance(particle1,particle2);
	b = (x1-x2)*(u1-u2)+(y1-y2)*(v1-v2);
	hoge = b*b-v_relative*v_relative*(r_relative*r_relative-pow(particle1.a+particle2.a,2.0));
	if(hoge > 0.0){
		tau += -(b+sqrt(hoge))/(v_relative*v_relative);
	}else{
		tau += T;
	}
	return tau;
}

double T_DWC(struct PARTICLE particle,double t,int j){//粒子と壁の衝突の時間を計算
	double tau = t;
	if(j==-1){//collision with RIGHT wall(-1)
		if(particle.u>0.0){
			tau += (Xmax-particle.a-particle.x)/particle.u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-2){//collision with LEFT wall(-2)
		if(particle.u<0.0){
			tau += (Xmin+particle.a-particle.x)/particle.u;
		}else{
			tau += 2.0*T;
		}
	}else if(j==-3){//collision with BOTTOM wall(-3)
		tau += (particle.v+sqrt(particle.v*particle.v+2*g*(particle.y-particle.a)))/g;
	}else if(j==-4){
		if(particle.x > particle.a){
			tau += -(particle.x-particle.a)/particle.u;
		}else if(particle.x < -particle.a){
			tau += (-particle.x-particle.a)/particle.u;
		}
	}else{//j=-5
		if(fabs(particle.x-X0) <= particle.a){
			tau += SolveQuarticEquation(g*g/4.0,-particle.v*g,-(particle.y-h)*g+particle.u*particle.u+particle.v*particle.v,2.0*((particle.x-X0)*particle.u+(particle.y-h)*particle.v),pow(particle.x-X0,2.0)+pow(particle.y-h,2.0)-particle.a*particle.a);
		}
	}
	if(tau < t){
		return 2.0*T;
	}else{
		return tau;
	}
}
double NextEvent(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,int j_current){//i_currentとj_currentのイベントを実際に行い，その時刻を返す関数
	double t = particle[i_current].event.time;
	Free_evolution(&particle[i_current],t-particle[i_current].tau);//i_currentの時間発展
	if(j_current >= 0){//Disk Disk Collision
		Free_evolution(&particle[j_current],t-particle[j_current].tau);//j_currentの時間発展
		G2(&particle[i_current],&particle[j_current]);//粒子同士の衝突処理
		if(r_distance(particle[i_current],particle[j_current]) < particle[i_current].a+particle[j_current].a-epsilon){
			printf("%d %d is too close!!\n",i_current,j_current);
			printf("distance = %lf\n",r_distance(particle[i_current],particle[j_current]));
		}
	}
	if(j_current < 0){//Disk Wall Collision
		G1(&particle[i_current],j_current);//壁との衝突処理
	}
	particle[i_current].event = Predictions(particle,cell,t,i_current);//i_currentのイベント更新
	CBT_update(node,particle[i_current].event.time,i_current);//i_currentのnodeアップデート
	if(j_current >= 0){//j_currentについても同様
		particle[j_current].event = Predictions(particle,cell,t,j_current);
		CBT_update(node,particle[j_current].event.time,j_current);
	}
	return t;
}

double t_cell_update(struct PARTICLE particle,int j_current,double t_cell_old,double *v_max){//セルの更新時刻を計算
	double t_cell,dt_cell;
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	if(j_current == -3){
		if(*v_max*(*v_max) < pow(particle.u,2.0)+pow(particle.v,2.0)){
			*v_max = sqrt(pow(particle.u,2.0)+pow(particle.v,2.0));
		}
	}
	dt_cell = (cell_length_y-2.0*a_large)/(2.0*(*v_max));
	t_cell = t_cell_old+dt_cell;
	return t_cell;
}

double Vmax(struct PARTICLE particle[N]){//最大速度を計算
	double v_max = 0.0;
	for(int i=0;i<N;i++){
		if(v_max*v_max < particle[i].u*particle[i].u+particle[i].v*particle[i].v){
			v_max = sqrt(particle[i].u*particle[i].u+particle[i].v*particle[i].v);
		}
	}
	return v_max;
}

void MaskUpdate(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],int i_current,double t){//同じマスクに含まれる粒子のイベントを更新
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	int cell_x = getcell_x(particle[i_current].x,cell_length_x) , cell_y = getcell_y(particle[i_current].y,cell_length_y),j;
	for(int c1=-1;c1<=1;c1++){
		for(int c2=-1;c2<=1;c2++){
			if((((cell_x+c1 >= 0) && (cell_x+c1 < N_cell_x)) && (cell_y+c2 >= 0)) && (cell_y+c2 < N_cell_y)){
				j = cell[cell_x+c1][cell_y+c2].first;
				while(j >= 0){
					if(particle[j].event.number_particle == i_current){
						particle[j].event = Predictions(particle,cell,t,j);
						CBT_update(node,particle[j].event.time,j);
					}
					j = particle[j].next;
				}
			}
		}
	}
}

double EEPGM(struct PARTICLE particle[N],struct CELL cell[N_cell_x][N_cell_y],struct NODE *node[n+1][2*p+2*q],double t,double *v_max){//すべての粒子に関して時間発展させたのちセルに登録し直す
	double cell_length_x = (Xmax-Xmin)/(double)N_cell_x,cell_length_y = (Ymax-Ymin)/(double)N_cell_y;
	double dt_cell,t_cell;
	double t_free,t_register,t_prediction,t_build,t_vmax;
	clock_t  start,end;
	start = clock();
	for(int i=0;i<N;i++){//現在時刻まで時間発展
		Free_evolution(&particle[i],t-particle[i].tau);
	}
	end = clock();
	t_free = (double)(end-start);
	
	start = clock();
	cell_register(particle,cell);//全粒子をセルに登録し直す
	end = clock();
	t_register = (double)(end-start);
	
	start = clock();
	for(int i=0;i<N;i++){//全粒子についてeventを計算し直す
		particle[i].event = Predictions(particle,cell,t,i);
	}
	end = clock();
	t_prediction = (double)(end-start);
	
	start = clock();
	CBT_build(node,particle);//CBTも最初から構成
	end = clock();
	t_build = (double)(end-start);
	
	start = clock();
	*v_max = Vmax(particle);
	end = clock();
	t_vmax = (double)(end-start);
	//printf("free:%lf,register:%lf,prediction:%lf,build:%lf,vmax:%lf\n",t_free,t_register,t_prediction,t_build,t_vmax);
	dt_cell = (cell_length_y-2.0*a_large)/(2.0*(*v_max));
	t_cell = t+dt_cell;
	return t_cell;
}

void SolveCubicEquation(double complex x[3],double a,double b,double c,double d){
	if(a == 0.0){
		printf("Error:a = 0.0\n");
		printf("This equation is NOT Cubic.\n");
	}
	double A,B,C,p_temp,q_temp,D;
	A = b/a;
	B = c/a;
	C = d/a;
	p_temp = B-A*A/3.0;
	q_temp = 2.0*A*A*A/27.0-A*B/3.0+C;
	D = q_temp*q_temp/4.0+p_temp*p_temp*p_temp/27.0;
	if(D < 0.0){//three real solutions
		double theta = atan2(sqrt(-D),-q_temp*0.5);
		x[0] = 2.0*sqrt(-p_temp/3.0)*cos(theta/3.0)-A/3.0;
		x[1] = 2.0*sqrt(-p_temp/3.0)*cos((theta+2.0*M_PI)/3.0)-A/3.0;
		x[2] = 2.0*sqrt(-p_temp/3.0)*cos((theta+4.0*M_PI)/3.0)-A/3.0;
	}else{//single real solution and two imaginary solutions(c.c)
		double u = Cuberoot(-q_temp*0.5+sqrt(D)),v = Cuberoot(-q_temp*0.5-sqrt(D));
		x[0] = u+v-A/3.0;
		x[1] = -0.5*(u+v)+sqrt(3.0)*0.5*1i*(u-v)-A/3.0;
		x[2] = -0.5*(u+v)-sqrt(3.0)*0.5*1i*(u-v)-A/3.0;
	}
}



double SolveQuarticEquation(double a,double b,double c,double d,double e){
	if(a == 0.0){
		printf("Error:a = 0.0\n");
		printf("This equation is NOT Quartic.\n");
	}
	
	double A,B,C,D,p_temp,q_temp,r_temp;
	A = b/a;
	B = c/a;
	C = d/a;
	D = e/a;
	p_temp = -6.0*pow(A/4.0,2.0)+B;
	q_temp = 8.0*pow(A/4.0,3.0)-2.0*B*A/4.0+C;
	r_temp = -3.0*pow(A/4.0,4.0)+B*pow(A/4.0,2.0)-C*A/4.0+D;
	double complex t_temp[3];
	SolveCubicEquation(t_temp,1.0,-p_temp,-4.0*r_temp,4.0*p_temp*r_temp-q_temp*q_temp);
	double t = creal(t_temp[0]);
	double complex m_temp = Squreroot(t-p_temp);//mは誤差なく実数or純虚数
	double tau = 2.0*T;
	double complex x[4];
	x[0] = (-m_temp+Squreroot(-t-p_temp+2.0*q_temp/m_temp))*0.5-A/4.0;
	x[1] = (-m_temp-Squreroot(-t-p_temp+2.0*q_temp/m_temp))*0.5-A/4.0;
	x[2] = (m_temp+Squreroot(-t-p_temp-2.0*q_temp/m_temp))*0.5-A/4.0;
	x[3] = (m_temp-Squreroot(-t-p_temp-2.0*q_temp/m_temp))*0.5-A/4.0;
	for(int i=0;i<4;i++){
		if(cimag(x[i]) == 0.0){//解が実
			if((creal(x[i]) < tau) && (creal(x[i]) > 0.0)){
				tau = creal(x[i]);
			}
		}
	}
	if(fabs(a*pow(tau,4.0)+b*pow(tau,3.0)+c*pow(tau,2.0)+d*tau+e) > epsilon){
		tau = 2.0*T;
	}
	return tau;
}



double Cuberoot(double x){
	if(x > 0.0){
		return pow(x,1.0/3.0);
	}else{
		return -pow(-x,1.0/3.0);
	}
}


double complex Squreroot(double complex x){
	double complex y;
	double r = sqrt(creal(x)*creal(x)+cimag(x)*cimag(x)),theta = atan2(cimag(x),creal(x));
	if(cimag(x) == 0.0){//t-pは実数だからかならずこちらで処理される
		if(creal(x) > 0.0){
			y = sqrt(r);
		}else{
			y = sqrt(r)*1i;
		}
	}else{
		if(theta < 0.0){
			theta += 2.0*M_PI;
		}
		double complex y = sqrt(r)*(cos(theta*0.5)+1i*sin(theta*0.5));
	}
	return y;
}


double m_cal(struct PARTICLE particle[N]){
	int N_left = 0,N_right;
	for(int i=0;i<N;i++){
		if(particle[i].x < X0){
			N_left += 1;
		}
	}
	N_right = N-N_left;
	double m = (intmax(N_left,N_right)-0.5*N)/(double)N;
	return m;
}

double intmax(int x,int y){
	if(x > y){
		return (double)x;
	}else{
		return (double)y;
	}
}


double bias(struct PARTICLE particle[N],char c[10]){
	double temp=0.0;
	if(strcmp(c,"small")==0){//small
		for(int i=0;i<N;i++){
			if(( strcmp(particle[i].size,"small")==0) && (particle[i].x < X0)){
				temp += 1.0;
			}
		}
	}else{//large
		for(int i=0;i<N;i++){
			if(( strcmp(particle[i].size,"large")==0) && (particle[i].x < X0)){
				temp += 1.0;
			}
		}
	}
	return temp;
}
