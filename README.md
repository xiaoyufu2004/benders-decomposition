# 经典Benders分解
## 1. 算法原理
**Benders分解**是行生成方法即**割平面法**的一种。割平面法的思想是，确定最优解可能不需要全部约束，因此可以逐行增加约束（割平面），减小可行域。

经典Benders分解针对**混合整数线性规划（MILP）**问题。对于一个MILP原问题，将其分解为**主问题**和**子问题**，并将原问题的变量分给这两个问题求解，一般将复杂变量或整数变量放在主问题中，将连续变量放在子问题中。先求解主问题，得到原问题下界，并将解传递给子问题。之后求解**对偶子问题**，若对偶子问题**有最优解**，则也得到子问题最优解，更新原问题上界，生成最优性割。若对偶子问题**无界**，则子问题无解，生成可行性割。将最优性割或可行性割加到主问题中，求解主问题并更新下界，迭代直到上下界很小为止。

## 2. 问题拆解
### 2.1 原问题：

原问题是一个最小化MILP问题， $\mathbf{x}$ 为整数变量， $\mathbf{y}$ 为连续变量， $\mathbf{x}$ 、 $\mathbf{y}$ 均为向量。

$$\min \quad \mathbf{c}^\top \mathbf{x} + \mathbf{d}^\top \mathbf{y} $$
$$\text{s.t.} \quad \mathbf{A} \mathbf{x} + \mathbf{B} \mathbf{y} \ge \mathbf{b} $$
$$ \mathbf{y} \ge 0 $$
$$ \mathbf{x} \in \mathbb{X}$$

### 2.2 主问题：

主问题的变量是原问题的**整数变量** $\mathbf{x}$ ，目标函数与原问题相同。主问题不包含原问题的约束条件，因此是**松弛主问题**。求解松弛主问题，得到原问题的一个下界。主问题中， $\theta$ 是辅助变量，用来“占位”表示“对子问题目标函数 $\phi(x)$ 的上界估计”，通过不断加割把 $\theta$ 收紧，最终逼得 $\theta=\phi(x)$ 。

$$\min \quad \mathbf{c}^\top \mathbf{x} +  \theta $$
$$\text{s.t.} \quad  \mathbf{x} \in \mathbb{X}$$


### 2.3 子问题：

将松弛主问题的求解结果带入子问题。子问题的变量是原问题的**连续变量** $\mathbf{y}$ ，目标函数 $\phi(\mathbf{x})=\mathbf{d}^\top \mathbf{y}$ 。

$$\phi(\mathbf{x}) = \min \quad  \mathbf{d}^\top \mathbf{y} $$
$$\text{s.t.} \quad  \mathbf{B} \mathbf{y} \ge \mathbf{b} - \mathbf{A} \bar{\mathbf{x}} $$
$$ \mathbf{y} \ge 0$$


### 2.4 对偶子问题：
在经典Benders分解中，取对偶的前提是子问题是一个**LP问题**。LP问题具有**强对偶性**，若对偶子问题有最优解，则对偶子问题最优解也必定是子问题最优解。

子问题目标函数是最小化，因此对偶问题目标函数是最大化：

$$\max \quad  (\mathbf{b} - \mathbf{A} \bar{\mathbf{x}})^\top \mathbf{u}$$
$$\text{s.t.} \quad  \mathbf{B}^\top \mathbf{u} \le \mathbf{d}$$
$$\mathbf{u} \ge 0$$

## 3. 可行性割与最优性割：
### 3.1 可行性割：

我们先引入**极射线**的概念。若对偶子问题**无界**（子问题无解），则其可行域存在某个方向（称为极射线），在这个方向上目标函数可以无限增加。极射线本质上是一种“逃离有界区域的方向”，描述了约束系统失效的方向性证据。我们利用这个极射线，构造一个**可行性割（feasibility cut）**，这个割将被引入主问题，以排除那些会导致子问题不可行的 $\mathbf{x}$ 值：

$$\mathbf{r}_t^\top (\mathbf{b} - \mathbf{A}\mathbf{x}) \le 0 \quad \text{for } t = 1, 2, \dots, T$$

### 3.2 最优性割

若对偶子问题在当前主问题解 $\bar{\mathbf{x}}$ 下有最优解，则说明子问题在该 $\bar{\mathbf{x}}$ 下是可行并有界的。此时，我们可以从中获得子问题的一个最优解 $\mathbf{y}^\ast$ ，从而构造出原问题的一个可行解 $(\bar{\mathbf{x}} , \mathbf{y}^\ast)$ ，其为原问题提供了一个**上界**，可用于更新当前原问题的最优值。

继续利用对偶理论，我们还可以从对偶子问题的最优对偶解 $\mathbf{u}^\ast$ 中提取关于子问题最优值函数 $\phi(x)$ 的**线性下界**：

$$
\phi(x) \ge (\mathbf{b} - \mathbf{A} \mathbf{x})^\top \mathbf{u}^\ast
$$

为在主问题中逼近 $\phi(x)$ ，我们令主问题中的辅助变量 $\theta$ 代表对子问题最优值的估计，并引入如下**最优性割（optimality cut）**，保证 $\theta$ 永远不小于子问题的真实最优值，从而保持目标函数的正确性（这里推导比较复杂）：

$$
\theta \ge (\mathbf{b} - \mathbf{A} \mathbf{x})^\top \mathbf{u}^\ast
$$

## 4 迭代流程

### 4.1 初始化
* 选取一个可行或启发式起点 $x^{(0)}$  
* 设上下界  

$$   \mathrm{UB}\leftarrow+\infty,\qquad \mathrm{LB}\leftarrow-\infty,\qquad k\leftarrow0   $$


### 4.2 循环迭代

#### 4.2.1 子问题  

$$\phi\!(x^{(k)})=\min_{y\ge0} \{d^{\top}y \bigm|By \ge b-Ax^{(k)}\}$$

* **若不可行**：取无界射线 $r_t$ ，生成可行性割

$$r_t^{\top}(b-Ax)\le 0$$

* **若可行且有界**：得对偶最优解 $u_s^{\star}$ ，生成最优性割，并更新上界  

$$  \theta \;\ge\; (b-Ax)^{\top}u_s^{\star}  $$

$$  \mathrm{UB}\;\leftarrow\;  \min\!\{\,\mathrm{UB},\;c^{\top}x^{(k)}+\phi\!(x^{(k)})\}  $$

#### 4.2.2 主问题（累积所有割）  

$$\min_{x\in\mathbb X,\;\theta}\; c^{\top}x + \theta \text{s.t.已生成的全部可行性割与最优性割}$$

得到新解 $(x^{(k+1)},\theta^{(k+1)})$ ，并更新下界  

$$\mathrm{LB}\;\leftarrow\;c^{\top}x^{(k+1)}+\theta^{(k+1)},\qquad
k\;\leftarrow\;k+1$$


### 4.3 终止
当 $\mathrm{UB}-\mathrm{LB}\le\varepsilon$ 时，  

$$x^{(k)},\;y^{\star}\!(x^{(k)})$$

即为原问题的近似最优解。
