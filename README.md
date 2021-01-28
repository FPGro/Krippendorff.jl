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
