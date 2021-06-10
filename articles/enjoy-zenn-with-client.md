---
title: "えんじょい"
emoji: "✨"
type: "tech" # tech: 技術記事 / idea: アイデア
topics: [FEM]
published: false
---
# 有限要素法
## 0.はじめに
有限要素法(FEM)をpythonで実装しました．本記事ではFEMの説明を行うとともにコードの説明を行っていきます．
私は現在大学生で来年には就活を控えています．「学生時代にうちこんだことはありますか？」と聞かれたときに「プログラミングを独学でしていました！」と答えるともりです．しかし面接官には，私がfor文を理解してドやっているのか，それともシステム管理やアプリ開発などおこなえるレベルなのかわかりません．よって今現在の私のレベルを伝えるべく執筆しました．

## 1.有限要素法とは
有限要素法とは，あらゆる形状の物体や荷重条件でも，応力やひずみなどの分布を計算できる数値解析手法の一つである．


## 2.仮想仕事の原理
一般に自然現象はエネルギーが最も低い状態になろうとする．この原理は材料に力や変位が加えられたときに生じる弾性変形にも適用できる．すなわち材料に力や変位を与えたときその材料のポテンシャルエネルギーが最小になる形に変形することが，材料にとって最も安定した状態であり，力学的に正しい解である．これはポテンシャルエネルギーの極小解を求めることを意味する．有限要素法では仮想仕事の原理を用いて各節点の変位を求めていく．

仮想仕事の原理を導出するためにまず，ポテンシャルエネルギーを求める，材料内に加わる表面力を$t$，物体力を$b$，変位を$u$とおく．

外力によるポテンシャルエネルギー$T_{\text{out}}$と材料内に蓄えられるポテンシャルエネルギー$T_{\text{in}}$の和が全ポテンシャルエネルギー$\Pi$である．$\Pi$を式で表すと

$$\Pi = T_{\text{in}} + T_{\text{out}}$$

$$= \frac{1}{2}\int_{\Gamma}^{\ }{(\sigma_{x}\varepsilon_{x} + \sigma_{y}\varepsilon_{y} + \sigma_{z}\varepsilon_{z} +}\tau_{\text{xy}}\gamma_{\text{xy}} + \tau_{\text{yz}}\gamma_{\text{yz}} + \tau_{\text{xz}}\gamma_{\text{xz}})dV$$

$$- \int_{S}^{\ }\left( ut_{x} + ut_{y} + ut_{z} \right)\text{dS}$$

$$\begin{matrix}
 - \int_{\Gamma}^{\ }{(ub_{x} + ub_{y} +}ub_{z})dV\\
\end{matrix}$$

$\Pi$の変化量$\delta\Pi = \delta T_{\text{in}} + \delta T_{\text{out}}$が0となるときの変位$u$の分布が，力がつり合っているときの解である．この状態を仮想仕事の原理とよび，式で表すと

$$\begin{matrix}
\int_{\Gamma}^{\ }{\left\{ \sigma \right\}^{T}\left\{ \text{δε} \right\} dV = \int_{S}^{\ }{\left\{ t \right\}^{T}\left\{ \text{δu} \right\}\text{dS}} + \int_{\Gamma}^{\ }{\left\{ b \right\}^{T}\left\{ \text{δu} \right\}\text{dV}}}\ \\
\end{matrix}$$

である．

## 3.要素と節点
2次元の有限要素法では材料を三角形または四角形要素に分割しそれぞれの要素の変位を求めていく． 本論文のモデルは4節点四角形要素で解析を行った．四角形要素は三角形要素に比べて同じ節点数でも要素数が少なくて済み，また要素内の応力分布，ひずみ分布を仮定しているので，滑らかな分布が得られる特徴があるからである．

## 4.全体剛性方程式
### 4.1.概要
2次元有限要素法では分割した要素に対して剛性方程式を求め，すべての要素についての剛性方程式を合わせることで全体剛性方程式を作り，各節点の変位を求める．得られた節点変位からひずみ，応力が求められる．4.2，4.3，4.4で剛性方程式を求めるために必要な知識を示し，4.5で全体剛性方程式を示す． 

### 4.2.*N*マトリックス
*N*マトリックスは形状関数と呼ばれ，要素内の変位$u$を節点変位$u_{m}$で表すときに用いられるものである．

要素内の変位*u*，*m*番目の要素の節点変位ベクトル*u~m~*，*m*番目の要素の*N*マトリックス*N~m~*の関係は

$$\left\{ \mathbf{u} \right\}\mathbf{=}\begin{pmatrix}
\mathbf{u} \\
\mathbf{v} \\
\end{pmatrix}\mathbf{=}\begin{bmatrix}
\mathbf{N}_{\mathbf{1}} & \mathbf{0} & \mathbf{N}_{\mathbf{2}} & \mathbf{0} & \mathbf{N}_{\mathbf{3}} & \mathbf{0} & \mathbf{N}_{\mathbf{4}} & \mathbf{0} \\
\mathbf{0} & \mathbf{N}_{\mathbf{1}} & \mathbf{0} & \mathbf{N}_{\mathbf{2}} & \mathbf{0} & \mathbf{N}_{\mathbf{3}} & \mathbf{0} & \mathbf{N}_{\mathbf{4}} \\
\end{bmatrix}\begin{pmatrix}
\mathbf{u}_{\mathbf{1}} \\
\mathbf{v}_{\mathbf{1}} \\
\mathbf{u}_{\mathbf{2}} \\
\mathbf{v}_{\mathbf{2}} \\
\mathbf{u}_{\mathbf{3}} \\
\mathbf{v}_{\mathbf{3}} \\
\mathbf{u}_{\mathbf{4}} \\
\mathbf{v}_{\mathbf{4}} \\
\end{pmatrix}$$

$$\begin{matrix}
\mathbf{=}\left\lbrack \mathbf{N}_{\mathbf{m}}\mathbf{(}\xi,\eta\mathbf{)} \right\rbrack\left\{ \mathbf{u}_{\mathbf{m}} \right\}\mathbf{\#}(2.3) \\
\end{matrix}$$

と表される．ここで*N~1~*，*N~2~*，*N~3~*，*N~4~*は次式で表される．

$$\begin{matrix}
N_{1} = \frac{1}{4}(1 - \xi)(1 - \eta)\#(2.4) \\
\end{matrix}$$

$$\begin{matrix}
N_{2} = \frac{1}{4}(1 + \xi)(1 - \eta)\#(2.5) \\
\end{matrix}$$

$$\begin{matrix}
N_{3} = \frac{1}{4}(1 + \xi)(1 + \eta)\#(2.6) \\
\end{matrix}$$

$$\begin{matrix}
N_{4} = \frac{1}{4}(1 - \xi)(1 + \eta)\#(2.7) \\
\end{matrix}$$

$\xi$および$\eta$について説明する．四角形要素は任意の形状をとりうるので$x$-$y$座標のまま計算を行うのは何かと不便な点が多いそこで図2.1(a)のように局所座標系$O$-$\text{ξη}$をとる．$O$-$\text{ξη}$は各辺の中点を通りその切片が1となるように定義する．子の座標系による要素は2.1(b)のように一片の長さが2の正方形に変換される．

### 4.3. *B*マトリックス

*B*マトリックスは要素内におけるひずみを節点変位で表すときに用いられるものである．

要素内のひずみε，節点変位*u~m~*，および*B*マトリックスの関係は

$$\left\{ \mathbf{\varepsilon} \right\}\mathbf{=}\begin{bmatrix}
\frac{\mathbf{\partial}\mathbf{N}_{\mathbf{1}}}{\mathbf{\partial x}} & \mathbf{0} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{2}}}{\mathbf{\partial x}} & \mathbf{0} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{3}}}{\mathbf{\partial x}} & \mathbf{0} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{4}}}{\mathbf{\partial x}} & \mathbf{0} \\
\mathbf{0} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{1}}}{\mathbf{\partial y}} & \mathbf{0} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{2}}}{\mathbf{\partial y}} & \mathbf{0} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{3}}}{\mathbf{\partial y}} & \mathbf{0} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{4}}}{\mathbf{\partial y}} \\
\frac{\mathbf{\partial}\mathbf{N}_{\mathbf{1}}}{\mathbf{\partial y}} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{1}}}{\mathbf{\partial x}} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{2}}}{\mathbf{\partial y}} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{2}}}{\mathbf{\partial x}} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{3}}}{\mathbf{\partial y}} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{3}}}{\mathbf{\partial x}} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{4}}}{\mathbf{\partial y}} & \frac{\mathbf{\partial}\mathbf{N}_{\mathbf{4}}}{\mathbf{\partial x}} \\
\end{bmatrix}\mathbf{\ }\begin{pmatrix}
\mathbf{u}_{\mathbf{1}} \\
\mathbf{v}_{\mathbf{1}} \\
\mathbf{u}_{\mathbf{2}} \\
\mathbf{v}_{\mathbf{2}} \\
\mathbf{u}_{\mathbf{3}} \\
\mathbf{v}_{\mathbf{3}} \\
\mathbf{u}_{\mathbf{4}} \\
\mathbf{v}_{\mathbf{4}} \\
\end{pmatrix}$$

$$\begin{matrix}
\mathbf{=}\mathbf{\lbrack B}_{\mathbf{m}}\mathbf{(}\xi,\eta\mathbf{)\rbrack}\left\{ \mathbf{u}_{\mathbf{m}} \right\}\mathbf{\#}(2.8) \\
\end{matrix}$$

と表すことができる．

### 4.4. 全体剛性行列

節点力{*F*}は表面力と物体力を各節点に加わる力として置き換えたものである．要素剛性マトリックス\[*K~m~*\]は*B*マトリックスおよび*D*マトリックスを用いて

$$\begin{matrix}
\left\lbrack K_{m} \right\rbrack = \int_{V_{m}}^{\ }{\left\lbrack B_{m} \right\rbrack^{T}\lbrack D\rbrack\left\lbrack B_{m} \right\rbrack\text{dV}}\#(2.10) \\
\end{matrix}$$

と表すことができる．全体剛性マトリックス\[*K*\]は\[*K~m~*\]の重ね合わせにより求めることができる．

　全体剛性方程式は，仮想仕事の原理を変形させた式である．全体剛性マトリックス\[*K*\]，変位{*U*}，節点力{*F*}を用いて

$$\begin{matrix}
\lbrack K\rbrack\left\{ U \right\} = \left\{ F \right\}\#(2.9) \\
\end{matrix}$$

と表される．


