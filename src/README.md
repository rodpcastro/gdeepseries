# Infinite-depth free-surface Green function

## Variables


## Green function, gradient and hessian matrix


## Series expansions for $F$, $F_X$ and $F_{XX}$

### Series expansion for $D_1$

$$
F = -2 e^{-Y} \mathrm{Ei}(Y) + 2 \sum_{n=1}^{\infty} \frac{1}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right],
$$

$$
\frac{\partial F}{\partial X} = \frac{4}{X} \sum_{n=1}^{\infty} \frac{n}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right],
$$

$$
\frac{\partial^2 F}{\partial X^2} = \frac{4}{X^2} \sum_{n=1}^{\infty} \frac{n(2n-1)}{n!^2} \left(-\frac{X^2}{4}\right)^n
\left[\sum_{m=1}^{2n} \frac{(m-1)!}{Y^m} - e^{-Y} \mathrm{Ei}(Y)\right].
$$

$\mathrm{Ei}$ is the exponential integral.

### Series expansion for $D_2$

$$
F = -\pi e^{-Y} Y_0(X) - \frac{2R}{X^2} +
\frac{2R}{X^2} e^{-Y} \sum_{n=0}^{\infty} \frac{(n+1)X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\,_2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right],
$$

$$
\frac{\partial F}{\partial X} = \pi e^{-Y} Y_1(X) + \frac{2Y}{XR} -
\frac{2R}{X} e^{-Y} \sum_{n=0}^{\infty} \frac{X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\,_2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right],
$$

$$
\frac{\partial^2 F}{\partial X^2} = \pi e^{-Y} \left[Y_0(X) - \frac{Y_1(X)}{X} \right] + 
\frac{2Y}{X^2 R} \left(\frac{Y^2}{R^2} - 2 + Y\right) -
\frac{2R}{X^2} e^{-Y} \sum_{n=1}^{\infty} \frac{n X^n}{n!}
\mathfrak{Re}\left[\mathrm{i}^{-n}\,_2F_1\left(\frac{1}{2}, -\frac{n}{2}; \frac{3}{2}; \frac{R^2}{X^2}\right)\right].
$$

$_2F_1$ is the Gauss hypergeometric function.

### Series expansion for $D_3$

$$
C_0 = 1 - e^{-Y}, \quad C_n = Y^{2n} - 2n Y^{2n-1} + 2n(2n-1)C_{n-1},
$$

$$
F = -\pi e^{-Y} [H_0(X)+Y_0(X)] - \frac{2(1-e^{-Y})}{X} -
\frac{2}{X} \sum_{n=1}^{\infty} \frac{(-1)^n (2n-1)!!}{2^n n! X^{2n}} C_n,
$$

$$
\frac{\partial F}{\partial X} = -2e^{-Y} +\pi e^{-Y} [H_1(X)+Y_1(X)] + \frac{2(1-e^{-Y})}{X^2} +
\frac{2}{X^2} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1) (2n-1)!!}{2^n n! X^{2n}} C_n,
$$

$$
\frac{\partial^2 F}{\partial X^2} = 
\pi e^{-Y} \left[H_0(X)+Y_0(X) - \frac{H_1(X)+Y_1(X)}{X}\right] - \frac{4(1-e^{-Y})}{X^3} -
\frac{2}{X^3} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1) (2n+2) (2n-1)!!}{2^n n! X^{2n}} C_n.
$$

$(2n-1)!!$ is the odd double factorial.

### Series expansion for $D_4$

$$
B_0 = \frac{1-e^{-Y}}{Y},
\quad B_1 = \left(\frac{1}{Y} - \frac{1}{Y^3}\right)e^{-Y} - 2\left(\frac{1}{Y^2} - \frac{1}{Y^3}\right),
\quad B_n = \frac{(-1)^{n+1}}{Y}e^{-Y} + \frac{2n(2n-1)}{Y^2}B_{n-1} + \frac{4n(n-1)}{Y^2}B_{n-2},
$$

$$
F = -\pi e^{-Y} [H_0(X)+Y_0(X)] - \frac{2(1-e^{-Y})}{R} -
\frac{2Y}{R} \sum_{n=1}^{\infty} \frac{(-1)^n (2n-1)!!}{2^n n!} \left(\frac{Y}{R}\right)^{2n} B_n,
$$

$$
\frac{\partial F}{\partial X} = -2e^{-Y} +\pi e^{-Y} [H_1(X)+Y_1(X)] + \frac{2X(1-e^{-Y})}{R^3} -
\frac{2XY}{R^3} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1)!!}{2^n n!} \left(\frac{Y}{R}\right)^{2n} B_n,
$$

$$
\frac{\partial^2 F}{\partial X^2} = \pi e^{-Y} \left[H_0(X)+Y_0(X) - \frac{H_1(X)+Y_1(X)}{X}\right] +
\frac{2(1-e^{-Y})}{R^3} \left(1 - \frac{3X^2}{R^2}\right) +
\frac{2Y}{R^3} \sum_{n=1}^{\infty} \frac{(-1)^n (2n+1)!!}{2^n n!}
\left(\frac{Y}{R}\right)^{2n} \left(1 - \frac{X^2 (2n+3)}{R^2}\right) B_n.
$$
