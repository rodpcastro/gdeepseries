# GDeepSeries
This Fortran library contains series expansions for the three-dimensional infinite-depth free-surface Green function, implemented according to the expressions defined in the work of *Shan* and *Wu* (2018).

An overview of the expressions implemented in this library can be found [here][src]. Alternatively, a Python implementation of the series expansions of $F$, $F_X$ and $F_{XX}$ can be found in this [blog post][rpcgds].

## Test
Simple tests were performed to compare the values computed by GDeepSeries against the numerical evaluation with Python libraries. The test results are summarized [here][test].

## Dependencies
GDeepSeries makes use of special functions implemented in [ColSpecF].

## References
1. Penghao Shan and Jiameng Wu. Highly precise approximation of free surface Green function and its high order derivatives based on refined subdomains. Brodogradnja, vol. 69, no. 1, pp. 53–70, 2018. <https://doi.org/10.21278/brod69104>

## License
The GDeepSeries code is distributed under the MIT License (see [LICENSE] file).

**Important Dependency Notice**: This project depends on [ColSpecF], which is distributed under multiple licenses. Review [ColSpecF's license][csf-license] for details.

<!-- links -->
[colspecf]: https://colspecf.readthedocs.io/
<<<<<<< HEAD
[src]: https://github.com/rodpcastro/gdeepseries/tree/main/src/README.md#infinite-depth-free-surface-green-function
=======
[src]: https://github.com/rodpcastro/gdeepseries/tree/main/src/README.md#expressions
>>>>>>> 7800dcfe729b9f92e341b0cd26172100663dc31b
[license]: https://github.com/rodpcastro/gdeepseries/blob/main/LICENSE
[csf-license]: https://github.com/rodpcastro/colspecf/blob/main/LICENSE
[rpcgds]: https://rodpcastro.github.io/posts/0013_3d_inf_depth_fsurface_gfunction/
[test]: https://github.com/rodpcastro/gdeepseries/blob/main/test/README.md#test-results
