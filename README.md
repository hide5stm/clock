# Granular Clock Simulator

粉体時計に関するシミュレーションを行いアニメーションと粒子の偏りの時間変化のグラフを生成する。

### How to simulate

```
$ make
gcc -Wall -Wno-unused-variable -pg -lm -O2 gif.c -o gif

$ ./gif
```

gif.cがsetting.txtを生成し、また"gnuplot gif.txt"を実行させることによりアニメーションが生成される。
たいていのパラメータの情報はsetting.txtに入れられるので自動で修正されるが、粒子の大きさやその比率はplot.txtで手作業で変更するしかない。

