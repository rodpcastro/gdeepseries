# Test results
Two simple tests were performed using [test-drive], a standard Fortran unit test framework. The first test consists of obtaining the absolute error of $F$, $F_X$ and $F_{XX}$, by comparing GDeepSeries' results and the integral expressions evaluated with [SciPy]. Below there are the maximum absolute errors for $F$ and its derivatives, followed by the absolute errors obtained for several points $(X, Y) \in [0.005, 40.0] × [0.005, 40.0]$.

$F$ maximum absolute error = $1.32 \times 10^{-9}$

$F_X$ maximum absolute error = $1.94 \times 10^{-9}$

$F_{XX}$ maximum absolute error = $6.42 \times 10^{-9}$

![error_f](https://raw.githubusercontent.com/rodpcastro/gdeepseries/refs/heads/main/test/error/errorf.svg)

![error_fx](https://raw.githubusercontent.com/rodpcastro/gdeepseries/refs/heads/main/test/error/errorfx.svg)

![error_fxx](https://raw.githubusercontent.com/rodpcastro/gdeepseries/refs/heads/main/test/error/errorfxx.svg)

A time comparison was also made to check how advantageous the series expansions are compared to the integral formulas. It was noted that the series expansion method implemented in Fortran is about 50 times faster than the integral expressions evaluated with SciPy.

The second test is the computation of the absolute error of the Green function, its gradient and hessian matrix, for four cases. The values obtained by GDeepSeries are compared against very precise approximations evaluated with [numdifftools]. The absolute error for each case also satistfies the condition of staying below $10^{-8}$.

> [!NOTE]
> Tests can be executed by running the command `fpm test` from the root directory.

## References
1. The Fortran Programming Language. 2024. test-drive: The simple testing framework. <https://github.com/fortran-lang/test-drive>
2. P. Virtanen et al. 2020. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nat. Methods 17, 3 (Mar.), 261–272. <https://doi.org/10.1038/s41592-019-0686-2>
3. P. A. Brodtkorb. 2025. numdifftools: Solve automatic numerical differentiation problems in one or more variables. <https://github.com/pbrod/numdifftools>.

<!-- links -->
[test-drive]: https://github.com/fortran-lang/test-drive/
[scipy]: https://scipy.org/
[numdifftools]: https://github.com/pbrod/numdifftools/
