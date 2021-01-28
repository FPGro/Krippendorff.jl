# TODO insert header
"""
I should totally document the module
"""
module Krippendorff

using Tables

const UNITSDEFAULT = :rows
const KNOWNMETRICS = [:nominal, :interval]

## common entry point

# warn that NamedTuples/Dicts of Vectors and Vectors of NamedTuples/Dicts 
# can satisfy the tables interface. See Krippendorff.istable
"""
    krippendorffs_alpha(input; units = $(UNITSDEFAULT), metric = :nominal, R = :discrete, silent = false)
    Krippendorff.alpha(input; kwargs)

Compute the Krippendorff's-α inter-rater reliability measure from the supplied input. The input will 
be checked to determine how to iterate over it, `missing` values will be handled automatically if present 
and the most efficient algorithm to compute α will be determined heuristically. By default, `Tables.jl` 
tables and table-like inputs are assumed to have columns representing raters and rows representing units.
See [`prepare_iterator`](@ref) for more information about the input requirements.

# Arguments
- `units::Union{Symbol,AbstractString}`: `rows` or `col(umn)s`, see [`prepare_iterator`](@ref) for explanation
- `metric`: a metric computing the squared distance between any pair of responses. Any of $(KNOWNMETRICS) or a custom function. Should satisfy `f(x,y) = [0 if x==y], [>0 otherwise]` but this is not enforced. See README for explanation of the default metrics.
- `R`: The space of possible responses. Either `:discrete` (relatively few possible responses, uses fast computation via coincidence matrix), `:continuous` (many possible responses up to continuous range, uses a slower but generically applicable algorithm with minimal allocation) or a precomputed value (implies discrete, but avoids searching for all unique values). If precomputed, it should be supplied as a tuple or vector of possible values (e.g. the output of unique(yourdata)), or as an appropriate range object where possible (slightly more efficient). 
- `silent::Bool`: set to disable all optional output (`@info` and `@warn`, doesn't affect error messages)

See also: [`compute_alpha_generical`](@ref), [`compute_alpha_with_coincidences`](@ref)
"""
function krippendorffs_alpha(input; units::Union{Symbol,AbstractString} = UNITSDEFAULT, metric = :nominal, R = :discrete, silent::Bool = false)
    units_iterator, elementtype, inputsize = _dispatch_units_iterator(input, Symbol(units); silent)
    implementation, possible_responses = _dispatch_responses(R, units_iterator, elementtype, inputsize)
    distance_metric_squared = _dispatch_metric(metric, elementtype; silent)
    return _compute_alpha(implementation, possible_responses, distance_metric_squared, units_iterator)
end

const alpha = krippendorffs_alpha # enables Krippendorff.alpha
export krippendorffs_alpha

## input type dispatch + find a meaningful elementtype + size of the whole thing

function _dispatch_units_iterator(input, units::Symbol; silent::Bool = false)
    units ∈ (:rows, :columns, :cols) || error("\"units\" must be one of :rows or :col(umn)s, got: $(units)")
    userows = units === :rows # use rows as units or not (use columns then)
    units, elementtype, sizeestimate = _units_eltype_and_size(input, userows; silent)
    # if missings are allowed, wrap every unit in a skipmissings
    prepared_iter, prepared_eltype = Missing<:elementtype ? ((skipmissing(u) for u in units), nonmissingtype(elementtype)) : (units, elementtype)
    return (prepared_iter, prepared_eltype, sizeestimate)
end

# use multiple dispatch to find an appropriate units_iterator, eltype and size for different types

function _units_eltype_and_size(input, userows::Bool; silent::Bool = false)
    # fallback, check Tables.jl interface via istable, else assume generic iterability 
    if Tables.istable(input)
        if userows
            rows = Tables.rows(input)
            schema = Tables.schema(rows)
            if isnothing(schema)
                # infer from first row
                if !silent
                    @info "Schema-less table detected, infering layout from first row only. This may lead to unexpected results."
                end
                firstrow, _ = iterate(Iterators.take(rows,1))
                elementtype = reduce((x,y)->Union{x,y}, (typeof(entry) for entry in firstrow), init=Base.Bottom)
                colnames = Tables.columnnames(firstrow)
                sizeestimate = length(rows)*length(colnames)
                let colnames = colnames
                    return (( [Tables.getcolumn(unit, name) for name in colnames] for unit in rows ), elementtype, sizeestimate)
                end
            else
                elementtype = reduce((x,y)->Union{x,y}, schema.types, init=Base.Bottom)        
                sizeestimate = length(rows)*length(schema.names)
                return (rows, elementtype, sizeestimate)
            end
        else # usecols
            cols = Tables.columns(input)
            schema = Tables.schema(cols)
            if isnothing(schema)
                # infer from first row
                if !silent
                    @info "Schema-less table detected, infering types from first column only. This may lead to unexpected results."
                end
                colnames = Tables.columnnames(columns)
                firstcol = Tables.getcolumn(columns, 1)
                elementtype = reduce((x,y)->Union{x,y}, (typeof(entry) for entry in firstcol), init=Base.Bottom)
                sizeestimate = length(firstcol)*length(colnames)
                let colnames = colnames
                    return ((Tables.getcolumn(columns, name) for name in colnames), elementtype, sizeestimate)
                end
            else
                elementtype = reduce((x,y)->Union{x,y}, schema.types, init=Base.Bottom)        
                sizeestimate = length(Tables.getcolumn(cols,1))*length(schema.names)
                let colnames = schema.names
                    return ((Tables.getcolumn(cols, name) for name in colnames), elementtype, sizeestimate)
                end
            end
        end 
    else
        # last fallback, wrap everything in values(...) to accomodate dicts
        units = (values(unit) for unit in values(input))
        sizeestimate = sum(length(unit) for unit in units)
        elementtype = reduce((x,y)->Union{x,y}, (eltype(unit) for unit in units), init=Base.Bottom)
        if !silent
            @info "Input has no notion of rows and columns, iterating it is expected to yield units."
        end
        return (units, elementtype, sizeestimate)
    end
end

function _units_eltype_and_size(input::AbstractMatrix, userows::Bool; silent::Bool = false)
    iter = userows ? eachrow(input) : eachcol(input)
    return (iter, eltype(input), *(length(axes(input,1)), length(axes(input,2))))
end

## dispatch possible responses

function _dispatch_responses(R, _,_,_)
    # if not a Symbol: assume R is precomputed and reasonable
    return (Coincidence(), R)
end

function _dispatch_responses(R::Union{Symbol,AbstractString}, itr, elementtype, inputsize)
    Rsym = Symbol(R)
    if Rsym === :continuous 
        return (Generic(), nothing)
    elseif Rsym === :discrete
        # get unique values of input, if there are lots of unique values 
        # in comparison to the input size, use generic version
        # also don't bother with coincidence etc if the input is tiny
        uniques = unique(Base.Iterators.flatten(itr))
        # TODO determine cutoffs
        if inputsize / length(uniques) > 0.5
            return (Generic(), nothing)
        else
            return (Coincidence(), optimize_possible_responses(uniques, elementtype))
        end
    else 
        error("Unknown input for R: $(R). Use either :discrete, :continuous or a precomputed set of responses.")
    end
end


# try to find a more efficient mapping strategy than dict(unique_values => indices)
# i.e. if the input looks like a range and is sufficiently dense, then make it a range
function _optimize_possible_reponses(uniquevalues, elementtype)
    # if extended to a range, how many empty indices may appear at worst
    allowedsparsity::Float64 = 0.5
    if elementtype <: Integer # relax to number?
        minR, maxR = extrema(uniquevalues)
        # TODO check for range with a constant non-one stepsize
        if length(uniquevalues)/(maxR-minR) >= allowedsparsity # assumes stepsize == 1
            if minR === oneunit(elementtype) # one-based is especially easy
                return Base.OneTo(maxR)
            else
                return UnitRange(minR, maxR)
            end
        else # making a range would introduce too many unused responses
            return uniquevalues
        end
    else # all non-numerics can not be directly mapped to indices anyways
        return uniquevalues
    end
end

# return a function that maps all possible values to indices into a coincidence matrix axis
_map_responses_to_indices(::Nothing) = nothing
_map_responses_to_indices(::Base.OneTo) = identity
_map_responses_to_indices(possible_responses::UnitRange) = Base.Fix2(-, possible_responses.start - 1)
function _map_responses_to_indices(possible_responses)
    let mapping = Dict(reverse.(collect(pairs(possible_responses))))
        return x -> mapping[x]
    end
end

## dispatch inputs for metric (and assert that it fits the elementtype)

function _dispatch_metric(metric::Union{Symbol, AbstractString}, elementtype; silent::Bool = false)
    m = Symbol(metric)
    if m === :nominal
        return (!=)
    elseif m === :interval
        elementtype <: Number || (!silent && @warn "Interval metric is unlikely to work with non-numeric responses.")
        return (x,y)->abs2(x-y)
    else 
        throw(ArgumentError("Metric $(m) is currently not known. Try to use one of: $(KNOWNMETRICS)"))
    end
end

function _dispatch_metric(metric, elt; _...) # takes everything that's not a symbol or string
    # assumes that it's a valid function defined on every pair of possible responses
    return metric
end

# from a distance metric and the set of possible responses mapping to the coincidence matrix indices, 
# compute a symmetric distance matrix once for later evaluation

function _compute_distance_matrix(possible_responses, mapping_function::F, distance_metric_squared) where {F}
    l = length(possible_responses)
    D = zeros(l,l)
    for rᵢ in possible_responses, rⱼ in possible_responses
        i = mapping_function(rᵢ)
        j = mapping_function(rⱼ)
        D[i,j] = distance_metric_squared(rᵢ, rⱼ)
    end
    return D
end

## computation backends

abstract type KrippendorffComputationStrategy end
struct Coincidence<:KrippendorffComputationStrategy end
struct Generic<:KrippendorffComputationStrategy end

# TODO docstrings
# using coincidence matrix
function _compute_alpha(::Coincidence, possible_responses, distance_metric_squared, units_iterator)
    mapping_function = _map_responses_to_indices(possible_responses)
    distances = _compute_distance_matrix(possible_responses, mapping_function, distance_metric_squared)
    l = length(possible_responses)
    C = zeros(l, l) # initialize coincidence matrix
    _fill_coincidence_matrix!(C, mapping_function, units_iterator)
    return _evaluate_coincidence_matrix(C, distances)
end

function _fill_coincidence_matrix!(C, mapping_function::F, units_iterator) where {F} # should force specialization
    for unit in units_iterator
        increment::Float64 = 1 / (count(x->true,unit) - 1) # count(...) because skipmissings have no length
        for i in eachindex(unit), j in eachindex(unit)
            if i<j # this assumes ordered indices, won't work for e.g. NamedTuples but I think this is rarely a problem
                cᵢ = mapping_function(unit[i])
                cⱼ = mapping_function(unit[j])
                C[ cᵢ, cⱼ ] += increment
                C[ cⱼ, cᵢ ] += increment
            end
        end
    end
end

function _evaluate_coincidence_matrix(C, distance_matrix)
    l = size(C,1)
    marginals = sum(C, dims=1)
    n = sum(marginals)
    disagreement_observed = sum( C[i,j] * distance_matrix[i,j] for i in 1:l for j in i+1:l)
    disagreement_expected = sum( marginals[i] * marginals[j] * distance_matrix[i,j] for i in 1:l for j in i+1:l)
    # it would be sufficient to evaluate just the upper or lower triangle of each, but I don't think this would be faster
    return 1. - (n-1)*(disagreement_observed / disagreement_expected)
end

# generic without coincidence matrix
function _compute_alpha(::Generic, _, distance_metric, units_iterator)
    disagreement_observed = 0.
    disagreement_expected = 0.
    n_all = 0
    # one huge pass over the data: always add pairs to the expected term and 
    # to the observed term only when within one unit
    # to avoid duplications, enumerate both units iterator and skip all redundant
    for (i_s, unit_s) in enumerate(units_iterator) 
        n_inunit = count(x->true,unit_s) # count(...) because skipmissings have no length
        n_inunit > 1 || continue
        n_all += n_inunit
        for (i_t, unit_t) in enumerate(units_iterator) 
            count(x->true,unit_t) > 1 || continue
            i_t < i_s && continue # skip backwards pairs of units
            if i_s == i_t
                # both iters in same unit, iter forward pairs only
                unitsum = sum( distance_metric(i,j) 
                                for (i_pos, i) in enumerate(unit_s) 
                                for (j_pos, j) in enumerate(unit_t) 
                                if i_pos ≥ j_pos) # count same position pairs once for consistency, should be equal to i_pos > j_pos unless the metric violates d(x,x)==0
                disagreement_observed += unitsum / (n_inunit-1)
                disagreement_expected += unitsum
            else
                # iters in different units, iter both full, add only to expected disagreement term
                disagreement_expected += sum( distance_metric(i,j) for i in unit_s for j in unit_t )
            end
        end
    end
    return 1. - ((n_all-1)*(disagreement_observed / disagreement_expected))
end

## utility stuff

# thin wrappers for backends with explanation 

"""
    compute_alpha_generical(units_iterator, squared_distance_metric)

Generic computation backend for Krippendorffs alpha, bypassing the creation of coincidence tables
entirely. The tradeof in this case is the necessity to iterate over all pairs of units. Thus, scaling
is typically much worse (in O(U²*R) for (U)nits and (R)aters) than when using coicidence matrices.
Nonetheless, this backend will be used by default if the number of possible responses is large 
compared to the size of the input, since a coincidence matrix can get huge in this case.
It is also preferable when all possible responses are not known beforehand or span a continuous
spectrum of values. 

# Arguments

- `units_iterator`: An iterable object which is assumed to yield units with no `missings`. If this is not the case, you can call [`prepare_iterator`](@ref) on it.
- `squared_distance_metric`: An object callable on any pair of possible responses in the supplied iterator. Should satisfy `f(x,y) = [0 if x==y], [>0 otherwise]` but this is not checked explicitly.
"""
compute_alpha_generical(units, metric) = _compute_alpha(Generic(), nothing, metric, units)

"""
    compute_alpha_with_coincidences(units, metric, possible_responses)

The default fast computation backend for Krippendorffs alpha. This will iterate over all units only once
and thus scales preferably if the number of differing possible responses is bounded. A coincidence matrix
is generated to keep track of the observed disagreement, while the expected disagreement will be computed
from the marginals of the coincidence matrix.
Be careful, if the number of possible responses is large, this backend may allocate a lot! 

# Arguments 

- `units_iterator`: An iterable object which is assumed to yield units with no `missings`. If this is not the case, you can call [`prepare_iterator`](@ref) on it.
- `squared_distance_metric`: An object callable on any pair of possible responses in the supplied iterator. Should satisfy `f(x,y) = [0 if x==y], [>0 otherwise]` but this is not checked explicitly.
- `R`: Space of possible responses. This is necessary to generate the coincidence matrix efficiently. It should be supplied as a tuple or vector of possible values (e.g. the output of unique(yourdata)), or as an appropriate range object where possible. 
"""
compute_alpha_with_coincidences(units, metric, R) = _compute_alpha(Coincidence(), R, metric, units)

"""
    prepare_iterator(input; units = $(UNITSDEFAULT))

Prepare an object for iteration by one of the `compute_alpha_...` functions. 
This involves determining how to iterate over units in the object and probing for the elementtype.
If appropriate, the `units` argument is used to determine the direction of iteration. This is however
not always possible. If no sense of direction is found, the heuristic will assume the input is already
a suitable iterator over units. Furthermore, if the input is found to contain [`missing`](@ref Base.missing) values
(or has them in it's eltype), all units will be wrapped in [`skipmissing`](@ref Base.skipmissing) automatically.

Since some seemingly unstructured iterables can satisfy the `Tables.jl` interface somewhat surprisingly
(Dict{Symbol,Vector} does, but not Dict{String,Vector} for example) and this may change the order of 
iteration implied, you can call the helper function [`Krippendorff.istable`](@ref) to see whether your
input looks like a `table` and if yes, how many rows and columns it appears to have.

# Arguments

- `input`: The input to be prepared. Should support generic iteration via [`iterate`](@ref Base.iterate), [`eachrow`](@ref Base.eachrow) or [`eachcol`](@ref Base.eachcol) or satisfy the `Tables.jl` interface as determined by [`Tables.istable`](@ref).
- `units::Union{Symbol,AbstractString}`: either `rows` or `col(umn)s`. This is used to determine how to iterate units in the input. For example, if the input iterator was a `Matrix`, `:rows` would make the function call `eachrow(input)` (and a little bit more). 
"""
function prepare_iterator(input, units::Union{Symbol,AbstractString} = UNITSDEFAULT) 
    units_iter, _, _ = _dispatch_units_iterator(input, Symbol(units))
    return units_iter
end

"""
    istable(input; IO = stdout)

A thin wrapper around [`Tables.istable`](@ref) that additionally prints how man rows and columns the input appears
to have when iterated through the `Tables.jl` interface. ([`Tables.columns`](@ref) specifically)
`IO` can be used to redirect the written output. Pass `IO=devnull` to supress output (making it equivalent to calling [`Tables.istable`](@ref))

# Examples

The `Tables.jl` interface assumes named columns and unnamed rows, which may lead to confusion 
if one wanted to pass a dictionary of rows for examples:

```jldoctest; setup = :(using Krippendorff, Tables)
julia> testmatrix = reshape(1:15, (3,5))
3×5 reshape(::UnitRange{Int64}, 3, 5) with eltype Int64:
 1  4  7  10  13
 2  5  8  11  14
 3  6  9  12  15

julia> Krippendorff.istable(testmatrix)
false

julia> Krippendorff.istable(Tables.table(testmatrix));
Input satisfies the Tables.jl table interface and appears to have 3 rows and 5 columns.

julia> testvectordict = Dict([k=>v for (k,v) in zip([:row1, :row2, :row3], eachrow(testmatrix))]); [println(entry) for entry in testvectordict];
:row1 => [1, 4, 7, 10, 13]
:row2 => [2, 5, 8, 11, 14]
:row3 => [3, 6, 9, 12, 15]

julia> Krippendorff.istable(testvectordict)
Input satisfies the Tables.jl table interface and appears to have 5 rows and 3 columns.
true
```
"""
function istable(input; IO = stdout)
    table = Tables.istable(input)
    if table
        itr = Tables.columns(input)
        cols = length(Tables.columnnames(itr))
        rows = cols==0 ? 0 : length(Tables.getcolumn(itr,1))
        println(IO, "Input satisfies the Tables.jl table interface and appears to have $(rows) rows and $(cols) columns.")
    end
    return table
end

end # module
