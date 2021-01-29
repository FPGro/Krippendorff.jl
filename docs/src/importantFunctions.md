# The Krippendorff Module

Currently, there are two computational backends: 
If the number of possible responses is bounded, one can construct a coincidence matrix with all observed pairs. Computing the observed disagreement is straightforward in this case, while the expected disagreement can easily be cpmputed from the marginals of the coincidence matrix. 
For all other cases, a generical computation strategy was implemented which avoids constructing
temporary tables entirely, but typically scales much worse. 

## Detailed API

```@docs
krippendorffs_alpha
compute_alpha_with_coincidences
compute_alpha_generical
prepare_iterator
istable
```