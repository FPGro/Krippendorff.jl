[![Build status](https://github.com/FPGro/Krippendorff/workflows/CI/badge.svg)](https://github.com/FPGro/Krippendorff/actions)
[![codecov.io](http://codecov.io/github/FPGro/Krippendorff/coverage.svg?branch=main)](http://codecov.io/github/FPGro/Krippendorff?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://FPGro.github.io/Krippendorff/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://FPGro.github.io/Krippendorff/dev)

# Krippendorff
Julia implementation of Krippendorff's Alpha Inter-rater reliability measure

## TODO:
- documentation
- check for correctness
- add remaining metrices
- add tests

## general plan
- common function entry point
    - multiple methods (Matrix, maybe Dict, Default to assuming a Tables.jl interface)
    - preprare units iterator for workhorse function
    - maybe keywords for customization
- workhorse function doing the actual computation
    - at least 2 versions
    - for "small" amount of possible values in the table, use fast method via contingency tables
    - direct computation fallback (slower)
    - liberal choice of metrices (nominal, ordinal, squared, interval, custom)
