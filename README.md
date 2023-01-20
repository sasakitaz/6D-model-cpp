# 6D-model-cpp
## これは何？
修士論文Part.Iに示した6自由度モデルで固有値と固有ベクトルを計算するプログラムです．
Hamiltonian行列を生成し，対角化を行います．場合によっては基底をmethaneの対称性A, E, Fに関して対称化してブロック対角化します．
Hamiltonian行列全体を生成・対角化するのではなく，量子数mに関して個別にブロック行列を生成するので，目的の量子数mを入力することに注意してください．

## 使い方
Eigen( https://eigen.tuxfamily.org/index.php?title=Main_Page )をダウンロードし，main.cppと同じ階層にEigenフォルダを入れてください．

parameter.hにポテンシャルパラメータ，その他質量などのパラメータを入力してください．
よい量子数mはここで入力します．
単位はcm^-1, Åです．

main.cppを実行すると対称化の方法を訊かれます．
ふつうの運用をする場合は0を入力してください．

### コンパイル環境
Windows10

g++ (MSYS2 MinGW x64)

LAPACK 3.10.1-1

Eigen 3.4.0

### コンパイルコマンド
g++ main.cpp -o main.exe -llapack -lblas -lm -O2

## ライセンス

## 内部実装

### main.cpp

### Parameter.h
ポテンシャルパラメータ，物理定数の宣言．
量子数mに依存する行列サイズの定義．

### Wigner3jSymbol.h
3-jシンボルの定義．

### WignerdSmall_pihalf.h
Symmetrize.hにおいて基底対称化に必要なd関数のリスト．
d関数の計算はothersフォルダのWignerdSmall.hで行った．

### tools.h
行列要素で頻繁に現れる道具を定義．
- 2つの3-jシンボルの積
- DVR

### MatrixElement.h
Hamiltonian行列要素を定義．

### Symmetrize.h
対称コマの回転子基底をA, E, Fで対称化．
対称化した基底でHamiltonian行列をブロック対角化して分割．

### diagonalization_dsyevd.h
実対称行列の対角化LAPACK呼び出し，基底対称化を行っていない場合に使用．

### diagonalization_zheevr.h
複素Hermite行列の対角化LAPACK呼び出し，基底対称化を行っている場合に使用．

### SymmetrizeDiagonalization.h
生成された行列をdiagonalization_dsyevd.h, diagonalization_zheevr.hを用いて対角化

### Eigenstate.h
固有値と固有ベクトルの表示，txtファイルへの書き出し．
