[![Build status](https://github.com/FPGro/Krippendorff.jl/workflows/CI/badge.svg)](https://github.com/FPGro/Krippendorff.jl/actions)
[![codecov.io](http://codecov.io/github/FPGro/Krippendorff.jl/coverage.svg?branch=main)](http://codecov.io/github/FPGro/Krippendorff.jl?branch=main)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://FPGro.github.io/Krippendorff.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://FPGro.github.io/Krippendorff.jl/dev)

# Krippendorff.jl

Julia implementation of Krippendorff's alpha Inter-rater reliability measure.

## Quickstart

```julia
using Krippendorff
data = ...
krippendorffs_alpha(data)
```

Should be sufficient in the majority of cases. Note that by default, Krippendorff assumes that rows 
in your data correspond to units and columns represent different raters. For more information, see 
the docs or consult the Julia help mode.

## TODOs:
- check for correctness: Basic checks have been done, but more and especially edge cases should be incorporated into the tests.
- add remaining metrices: At least ordinal, ratio, (bi)polar and circular metrics are still missing.
- add parallelism: The variant using a coincidence matrix should be relatively easy to parallelize. I'd like to play around with Transducers a little and may add a third fast backend eventually.
- find other irr implementations: A whole package for a single metric seems a little out of proportion. If more similar implementations exist (different metrices, etc.), it may be worthwile to collect them into a more wholesome package. (feel free to contact me in this direction)
- preparation step takes considerable time, make faster

Extra: take a look at the *About* page in the docs to learn why this package was created.