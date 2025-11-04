# GDeepSeries
This Fortran library contains series expansions for the three-dimensional infinite-depth free-surface Green function, implemented according to the expressions defined in the work of *Shan* and *Wu* (2018).

A quick overview of the series expansions and their implementation in Python can be found [here][rpcgds].

# Test
A simple test was performed to compare the values computed by GDeepSeries against the numerical evaluation of the integral expressions. The test results are summarized [here][test].

## Dependencies
GDeepSeries makes use of special functions implemented in [ColSpecF].

## References
1. Penghao Shan and Jiameng Wu. Highly precise approximation of free surface Green function and its high order derivatives based on refined subdomains. Brodogradnja, vol. 69, no. 1, pp. 53–70, 2018. <https://doi.org/10.21278/brod69104>

## License
The GDeepSeries code is distributed under the MIT License (see [LICENSE] file).

**Important Dependency Notice**: This project depends on [ColSpecF], which is distributed under multiple licenses. Review [ColSpecF's license][csf-license] for details.

<!-- links -->
[colspecf]: https://colspecf.readthedocs.io/
[license]: https://github.com/rodpcastro/gdeepseries/blob/main/LICENSE
[csf-license]: https://github.com/rodpcastro/colspecf/blob/main/LICENSE
[rpcgds]: https://rodpcastro.github.io/posts/0013_3d_inf_depth_fsurface_gfunction/
[test]: https://github.com/rodpcastro/gdeepseries/blob/main/test/README.md
