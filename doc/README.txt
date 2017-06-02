
KSHELL : nuclear shell model calculation code for K-computer era

KSHELL User's guide   2013/10/22
       東京大学原子核科学研究センター 清水 則孝 
       Noritaka SHIMIZU, shimizu@cns.s.u-tokyo.ac.jp

Thick restrat Lanczos法による M-scheme 殻模型対角化計算コードです。
OXBASHライクなシェルを備えており、
エネルギー準位、E2/M1/E1遷移確率を求めることができます。
Linux PC上で簡単に使え、OpenMPスレッド並列に対応しているので
CPUコアを有効に使えます。
同じ使い勝手で、MPI+OpenMPハイブリッド並列による大規模並列計算にも
使えることが特徴です。
Ni56 pf-shell KB3相互作用の基底状態(10億次元の行列対角化に相当)を、
145秒で計算できます。(＠東京大学FX10 512ノード、8192コア)



情報提供、バグ報告、改良、.sntファイルの寄贈など、清水まで
お寄せください。( shimizu@cns.s.u-tokyo.ac.jp ) 
このソフトウェアは無保証です。
いかなる形においてもこのソフトウェアを利用する際は、
下記の文献を引用してください。

N. Shimizu, ``Nuclear shell-model code for massive parallel computation, KSHELL'',
arXiv:1310.5431 [nucl-th] (2013).



--- 
必須動作環境 
* Fortran compiler (Fortran 95 + ISO-TR15581拡張)
* BLAS / LAPACK library
* python ver.2.4 or later ( ver.3系列は未対応 )

オプション
* OpenMP
* MPI-2 library (MPI-IO)
* C compiler (経過時間計測のみに使用)

動作確認済環境： 
 Intel Fortran ver.11.1 + Intel MKL
 GNU Fortran 4.6.3 + lapack
 富士通Fortran + SSLII on FX10/京計算機 )

新しめの Linux OS + Intel Fortran compiler を推奨します。
( Windows 7 + cygwin + gfortran でも一応動きましたが、
  OpenMPを有効にするとかえって遅くなりました。
  openmpi でコアの数だけプロセスを立てることをお勧めします。)

--- 
インストール・使い方

tar xvzf kshell.tar.gz 
cd kshell/src

Makefile を環境にあうよう編集(デフォルトは Intel Fortran)

make single # for single processor
または、
make mpi    # for MPI parallel

alias kshell_ui.py=`pwd`/../bin/kshell_ui.py

これでインストールは終了です。
その後、新しいディレクトリを作成、移動して、
kshell_ui.py を実行することで、
OXBASHライクな対話型インターフェースにより
シェルスクリプトを生成します。
生成したスクリプトを実行すると、計算が始まります。
計算終了後、summary_hoge.txt にエネルギー、
遷移確率の計算結果の要約が保存されます。
Q-momentなどのより詳細な情報は、
summary_hoge.txtに log-file として記されている
ファイル log.hoge を確認してください。

模型空間・有効相互作用は独自の .snt 形式に格納されています。
kshell/snt に、いくつかの相互作用が例として含まれていますし、
自作することも可能です。
また、nushell2snt.py を用いることにより、
Nushell(or OXBASH)用の模型空間ファイルと相互作用ファイルを
 .snt 形式に変換出来ます。たとえば、
kshell/bin/nushell2snt.py sd.sp usda.int
とすると、フォーマットが変換された sd_usda.snt が作成されます。



Thanks to 宇都野穣(原子力機構)、水崎高浩(専修大)、本間道雄(会津大)、
  角田佑介(東大理)、角田直文(東大理)      (敬称略)



以下、細かい解説
-------------------------------------------------------------------------
制限事項 1.
   Intel Fortran compilerなど integer(16)をサポートしていない
   Fortran compiler を用いる場合、
   陽子あるいは中性子の一粒子状態数の上限が62。
   gfortran を用いて kmbit = 16 でコンパイルすれば、上限は126。

制限事項 2.
   J-射影なしで M=0 状態を計算すると、同じJ状態間のM1遷移確率と、
   磁気モーメントを求めることができません。(表示されません)
   これは、Clebsch-Gordan係数の性質 <J, 0, 1, 0 | J, 0> = 0.0 によるものです。
   必要な場合、J-射影計算、あるいは、M=1状態を計算してください。

-------------------------------------------------------------------------
相互作用ファイル .snt の解説
  一粒子軌道の定義、1体相互作用行列要素、2体相互作用行列要素と続きます。
  "!"から始まる行はコメントとみなされます。

例：w.snt  
! model space
! 陽子軌道数 中性子軌道数 閉殻陽子数 閉殻中性子数
   3   3     8   8
! 軌道番号 n   l   j  tz
   1       0   2   3  -1    !  1 = p 0d_3/2
   2       0   2   5  -1    !  2 = p 0d_5/2
   3       1   0   1  -1    !  3 = p 1s_1/2
   4       0   2   3   1    !  4 = n 0d_3/2
   5       0   2   5   1    !  5 = n 0d_5/2
   6       1   0   1   1    !  6 = n 1s_1/2
! 1体相互作用
! 行数 method1
     6   0
!  i   j      <i|V|j>
   1   1      1.64658
   2   2     -3.94780
   3   3     -3.16354
   4   4      1.64658
   5   5     -3.94780
   6   6     -3.16354
! 2体相互作用
!  行数 method2  A   mass dependence factor
   158   1  18  -0.30000
!  i   j   k   l     J      <i,j| V | k,l>_J
   1   1   1   1     0     -2.18450
   1   1   1   1     2     -0.06650
... 158行続く 
TBMEは <i,j| V | k,l>_J で指定し、アイソスピンは仮定しない。
TBMEの質量依存性は、
  method2 = 0   質量依存性なし
  method2 = 1   (質量数/A)^(mass dependence factor)


-----------------------------------------------------------------------



相互作用ファイル：jisp_3shl_hw25.snt

no core 計算の場合、２体相互作用と運動エネルギー部分をわけてとりあつかう
! n_jorb(1), n_jorb(2), ncore, zcore
! index,  n,  l,  j, tz
   6   6     0   0
    1     0   0   1  -1
    2     0   1   1  -1
    3     0   1   3  -1
    4     1   0   1  -1
    5     0   2   3  -1
    6     0   2   5  -1
    7     0   0   1   1
    8     0   1   1   1
    9     0   1   3   1
   10     1   0   1   1
   11     0   2   3   1
   12     0   2   5   1
! interaction
! num, method=10,  hbar_omega
! method=10の時、運動エネルギーの行列要素であると仮定して
! hw*(A-1)/A の依存性ファクターをかける。
! エルミートを仮定して i<j
  i  j     <i|T_1b|j>
  14  10      25.00000000
  1   1      0.75000000
  1   4      0.61237244   ! 非対角要素 0s1/2-1s1/2
  2   2      1.25000000
  3   3      1.25000000
  4   4      1.75000000
  5   5      1.75000000
  6   6      1.75000000
  7   7      0.75000000
  7  10      0.61237244   ! 非対角要素
  8   8      1.25000000
  9   9      1.25000000
 10  10      1.75000000
 11  11      1.75000000
 12  12      1.75000000
! two-body int.
! num, method=10,  hbar_omega
! method=10の時、7番目のカラムを運動エネルギーの行列要素であると仮定して
! hw/A の依存性ファクターをかける。
! i   j   k   l   JJ     <ij|V_JISP|kl>_J  <ij|T_2b|kl>_J
        890  10      25.00000000
  1   1   1   1    0      -11.30897000      -0.00000000
  1   1   1   4    0       -4.61683600      -0.00000000
  1   1   2   2    0        2.66553100      -0.50000000
...



