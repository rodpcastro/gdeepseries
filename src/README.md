# Infinite-depth free-surface Green function
The infinite-depth free-surface Green function represents the spatial component of a velocity potential induced at a field point $p(x, y, z)$ by a pulsating source point $q(\xi, \eta, \zeta)$. Before presenting the Green function expressions, we define the following quantities

$$
\begin{split}
r = \sqrt{(x-\xi)^2 + (y-\eta)^2} &, \quad Z = z + \zeta, \\
X = k_0 r, \quad Y = -k_0 Z &, \quad \bar{Y} = Y + 2 k_0 \zeta, \\
\bar{R} = \sqrt{X^2 + \bar{Y}^2} &, \quad R = \sqrt{X^2 + Y^2}. \\
\end{split}
$$

where $k_0$ is the infinite-depth wave number. Between $p$ and $q$, $r$ is the horizontal distance and $\bar{R}$ is the Euclidian distance in dimensionless form. $|Z|$ is the vertical distance between $p$ and the image of $q$ with respect to the free-surface, while $R$ is the corresponding Euclidian distance in dimensionless form. $X$ and $Y$ are non-negative dimensionless cylindrical coordinates.

The infinite-depth free-surface Green function $G_\infty$, together with its first and second order derivatives, can be expressed as a funcion of $X$ and $Y$ as follows:

$$
G_\infty = \frac{k_0}{\bar{R}} + \frac{k_0}{R} + k_0 F + 2\mathrm{i}\pi s k_0 e^{-Y} J_0(X),
$$

$$
\frac{\partial G_\infty}{\partial X} = -\frac{k_0 X}{\bar{R}^3} - \frac{k_0 X}{R^3} + k_0 F_X - 2\mathrm{i}\pi s k_0 e^{-Y} J_1(X),
$$

$$
\frac{\partial G_\infty}{\partial Y} = -\frac{k_0 \bar{Y}}{\bar{R}^3} - \frac{k_0 Y}{R^3} + k_0 F_Y - 2\mathrm{i}\pi s k_0 e^{-Y} J_0(X),
$$

$$
\frac{\partial^2 G_\infty}{\partial X^2} = \frac{3 k_0 X^2}{\bar{R}^5} + \frac{3 k_0 X^2}{R^5} - \frac{k_0}{\bar{R}^3} - \frac{k_0}{R^3} + k_0 F_{XX} - 2\mathrm{i}\pi s k_0 e^{-Y} \left(J_0(X) - \frac{J_1(X)}{X}\right),
$$

$$
\frac{\partial^2 G_\infty}{\partial Y^2} = \frac{3 k_0 \bar{Y}^2}{\bar{R}^5} + \frac{3 k_0 Y^2}{R^5} - \frac{k_0}{\bar{R}^3} - \frac{k_0}{R^3} + k_0 F_{YY} + 2\mathrm{i}\pi s k_0 e^{-Y} J_0(X),
$$

$$
\frac{\partial^2 G_\infty}{\partial X \partial Y} = \frac{3 k_0 X \bar{Y}}{\bar{R}^5} + \frac{3 k_0 X Y}{R^5} + k_0 F_{XY} + 2\mathrm{i}\pi s k_0 e^{-Y} J_1(X),
$$

where $J_n$ is the $n$-th order Bessel function of the first kind. The variable $s=-1$ when the time component is $e^{\mathrm{i} \omega t}$, and $s=+1$ when the time component is $e^{-\mathrm{i} \omega t}$, where $\omega$ is the pulsating source frequency. The function $F$ and its derivatives are a function of $X$ and $Y$ given below

$$
F = -\pi e^{-Y} [H_0(X) + Y_0(X)] - 2 e^{-Y} \int_{0}^{Y} e^t (X^2+t^2)^{-\frac{1}{2}} \thinspace dt,
$$

$$
\frac{\partial F}{\partial X} = -2 e^{-Y} + \pi e^{-Y} [H_1(X) + Y_1(X)] +
2 X e^{-Y} \int_{0}^{Y} e^t (X^2+t^2)^{-\frac{3}{2}} \thinspace dt,
$$

$$
\frac{\partial F}{\partial Y} = -\frac{2}{R} - F,
$$

$$
\frac{\partial^2 F}{\partial X^2} = \frac{1}{3} e^{-Y} X + \frac{\pi}{2} e^{-Y} [H_0(X) + Y_0(X) - H_2(X) - Y_2(X)] +
2 e^{-Y} \int_{0}^{Y} e^t (X^2+t^2)^{-\frac{5}{2}} (t^2-2X^2) \thinspace dt,
$$

$$
\frac{\partial^2 F}{\partial Y^2} = \frac{2 Y}{R^3} + \frac{2}{R} + F,
$$

$$
\frac{\partial^2 F}{\partial X \partial Y} = \frac{2 X}{R^3} - F_X,
$$

where $H_n$ is the $n$-th order Struve function and $Y_n$ is the $n$-th order Bessel function of the second kind. From the equations above, it is noticeable that the derivatives with respect to $Y$ are related to $F$ and $F_X$, so the use of the series exansion method is focused only on $F$, $F_X$ and $F_{XX}$, as described in the following subtopic.

## Series expansions for $F$, $F_X$ and $F_{XX}$
*Shan & Wu (2018)* derived series expansions for $F$, $F_X$ and $F_{XX}$ for four different regions of the first quadrant of the $XY$ plane. The following image depicts these regions $D_i$, for $i=1,\ldots,4$ and the number of terms $N_i$ that are used to approximate $F$ and its derivatives in each region.

<!-- https://raw.githubusercontent.com/rodpcastro/gdeepseries/refs/heads/main/src/subdomains.svg -->

<p align="center">
<img src="https://raw.githubusercontent.com/rodpcastro/gdeepseries/refs/heads/gtest/src/subdomains.svg" alt="Subdomains">
</p>

### Series expansion in $D_1$

$$
F = -2 e^{-Y} \mathrm{Ei}(Y) + 2 \sum_{n=1}^{N_1} \frac{1}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right],
$$

$$
\frac{\partial F}{\partial X} = \frac{4}{X} \sum_{n=1}^{N_1} \frac{n}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right],
$$

$$
\frac{\partial^2 F}{\partial X^2} = \frac{4}{X^2} \sum_{n=1}^{N_1} \frac{n(2n-1)}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right].
$$

$\mathrm{Ei}$ is the exponential integral.

### Series expansion in $D_2$

$$
F = -\pi e^{-Y} Y_0(X) - \frac{2R}{X^2} +
\frac{2R}{X^2} e^{-Y} \sum_{n=0}^{N_2} \frac{(n+1)X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\thinspace _2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right],
$$

$$
\frac{\partial F}{\partial X} = \pi e^{-Y} Y_1(X) + \frac{2Y}{XR} -
\frac{2R}{X} e^{-Y} \sum_{n=0}^{N_2} \frac{X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\thinspace _2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right],
$$

$$
\frac{\partial^2 F}{\partial X^2} = \pi e^{-Y} \left[Y_0(X) - \frac{Y_1(X)}{X} \right] + 
\frac{2Y}{X^2 R} \left(\frac{Y^2}{R^2} - 2 + Y\right) -
\frac{2R}{X^2} e^{-Y} \sum_{n=1}^{N_2} \frac{n X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\thinspace _2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right].
$$

$_2F_1$ is the Gauss hypergeometric function.

### Series expansion in $D_3$

$$
C_0 = 1 - e^{-Y}, \quad C_n = Y^{2n} - 2n Y^{2n-1} + 2n(2n-1)C_{n-1},
$$

$$
F = -\pi e^{-Y} [H_0(X)+Y_0(X)] - \frac{2(1-e^{-Y})}{X} -
\frac{2}{X} \sum_{n=1}^{N_3} \frac{(-1)^n (2n-1)!!}{2^n n! X^{2n}} C_n,
$$

$$
\frac{\partial F}{\partial X} = -2e^{-Y} +\pi e^{-Y} [H_1(X)+Y_1(X)] + \frac{2(1-e^{-Y})}{X^2} +
\frac{2}{X^2} \sum_{n=1}^{N_3} \frac{(-1)^n (2n+1) (2n-1)!!}{2^n n! X^{2n}} C_n,
$$

$$
\frac{\partial^2 F}{\partial X^2} = 
\pi e^{-Y} \left[H_0(X)+Y_0(X) - \frac{H_1(X)+Y_1(X)}{X}\right] - \frac{4(1-e^{-Y})}{X^3} -
\frac{2}{X^3} \sum_{n=1}^{N_3} \frac{(-1)^n (2n+1) (2n+2) (2n-1)!!}{2^n n! X^{2n}} C_n.
$$

$(2n-1)!!$ is the odd double factorial.

### Series expansion in $D_4$

$$
B_0 = \frac{1-e^{-Y}}{Y},
\quad B_1 = \left(\frac{1}{Y} - \frac{1}{Y^3}\right)e^{-Y} - 2\left(\frac{1}{Y^2} - \frac{1}{Y^3}\right),
\quad B_n = \frac{(-1)^{n+1}}{Y}e^{-Y} + \frac{2n(2n-1)}{Y^2}B_{n-1} + \frac{4n(n-1)}{Y^2}B_{n-2},
$$

$$
F = -\pi e^{-Y} [H_0(X)+Y_0(X)] - \frac{2(1-e^{-Y})}{R} -
\frac{2Y}{R} \sum_{n=1}^{N_4} \frac{(-1)^n (2n-1)!!}{2^n n!} \left(\frac{Y}{R}\right)^{2n} B_n,
$$

$$
\frac{\partial F}{\partial X} = -2e^{-Y} +\pi e^{-Y} [H_1(X)+Y_1(X)] + \frac{2X(1-e^{-Y})}{R^3} -
\frac{2XY}{R^3} \sum_{n=1}^{N_4} \frac{(-1)^n (2n+1)!!}{2^n n!} \left(\frac{Y}{R}\right)^{2n} B_n,
$$

$$
\frac{\partial^2 F}{\partial X^2} = \pi e^{-Y} \left[H_0(X)+Y_0(X) - \frac{H_1(X)+Y_1(X)}{X}\right] +
\frac{2(1-e^{-Y})}{R^3} \left(1 - \frac{3X^2}{R^2}\right) +
\frac{2Y}{R^3} \sum_{n=1}^{N_4} \frac{(-1)^n (2n+1)!!}{2^n n!}
\left(\frac{Y}{R}\right)^{2n} \left(1 - \frac{X^2 (2n+3)}{R^2}\right) B_n.
$$

## Green function gradient and Hessian
The derivatives with respect to the field point coordinates are computed through chain rule as follows

$$
\frac{\partial G_\infty}{\partial x} = \frac{\partial G_\infty}{\partial X} \frac{\partial X}{\partial x},
$$

$$
\frac{\partial G_\infty}{\partial y} = \frac{\partial G_\infty}{\partial X} \frac{\partial X}{\partial y},
$$

$$
\frac{\partial G_\infty}{\partial z} = \frac{\partial G_\infty}{\partial Y} \frac{\partial Y}{\partial z},
$$

$$
\frac{\partial^2 G_\infty}{\partial x^2} = \frac{\partial G_\infty}{\partial X} \frac{\partial^2 X}{\partial x^2} + 
\frac{\partial^2 G_\infty}{\partial X^2} \left(\frac{\partial X}{ \partial x} \right)^2,
$$

$$
\frac{\partial^2 G_\infty}{\partial y^2} = \frac{\partial G_\infty}{\partial X} \frac{\partial^2 X}{\partial y^2} + 
\frac{\partial^2 G_\infty}{\partial X^2} \left(\frac{\partial X}{ \partial y} \right)^2,
$$

$$
\frac{\partial^2 G_\infty}{\partial z^2} = \frac{\partial G_\infty}{\partial Y} \frac{\partial^2 Y}{\partial z^2} + 
\frac{\partial^2 G_\infty}{\partial Y^2} \left(\frac{\partial Y}{ \partial z} \right)^2,
$$

$$
\frac{\partial^2 G_\infty}{\partial y \partial z} = 
\frac{\partial^2 G_\infty}{\partial X \partial Y} \frac{\partial X}{\partial y} \frac{\partial Y}{\partial z},
$$

$$
\frac{\partial^2 G_\infty}{\partial x \partial z} = 
\frac{\partial^2 G_\infty}{\partial X \partial Y} \frac{\partial X}{\partial x} \frac{\partial Y}{\partial z},
$$

$$
\frac{\partial^2 G_\infty}{\partial x^2} = \frac{\partial G_\infty}{\partial X} \frac{\partial^2 X}{\partial x \partial y} + 
\frac{\partial^2 G_\infty}{\partial X^2} \frac{\partial X}{ \partial x} \frac{\partial X}{ \partial y}.
$$

The derivatives of $X$ and $Y$ with respect to $x$, $y$ and $z$ are given below.

$$
\Delta u = x - \xi, \quad \Delta v = y - \eta, \quad r = \sqrt{{\Delta u}^2 + {\Delta v}^2}.
$$

$$
\frac{\partial X}{\partial x} = \frac{k_0 \Delta u}{r},
$$

$$
\frac{\partial X}{\partial y} = \frac{k_0 \Delta v}{r},
$$

$$
\frac{\partial^2 X}{\partial x^2} = \frac{k_0 {\Delta u}^2}{r^3},
$$

$$
\frac{\partial^2 X}{\partial y^2} = \frac{k_0 {\Delta v}^2}{r^3},
$$

$$
\frac{\partial^2 X}{\partial x \partial y} = -\frac{k_0 \Delta u \Delta v}{r^3},
$$

$$
\frac{\partial Y}{\partial z} = -k_0,
$$

$$
\frac{\partial^2 Y}{\partial z^2} = 0.
$$

## Reference
1. Penghao Shan and Jiameng Wu. Highly precise approximation of free surface Green function and its high order derivatives based on refined subdomains. Brodogradnja, vol. 69, no. 1, pp. 53–70, 2018. <https://doi.org/10.21278/brod69104>
