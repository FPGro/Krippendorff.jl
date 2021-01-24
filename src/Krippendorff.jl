module Krippendorff

using Tables: istable, rows, Columns

# TODO:
# - finish first draft
# - documentation
# - deal with missings (check eltype for Missing<:Eltype probably and add skipmissings)
# - check for correctness
# - add remaining metrices

# general plan
# - common function entry point
#     - multiple methods (Matrix, maybe Dict, Default to assuming a Tables.jl interface)
#     - preprare units iterator for workhorse function
#     - maybe keywords for customization
# - workhorse function doing the actual computation
#     - at least 2 versions
#     - for "small" amount of possible values in the table, use fast method via contingency tables
#     - direct computation fallback (slower)
#     - liberal choice of metrices (nominal, ordinal, squared, interval, custom)

## common entry point

# TODO docstrings
# input = input
# units = :rows, "rows", :columns, "columns" whether units are rows or columns 
# metric = one of (String or Symbol) nominal, ordinal, interval, ratio,...? or a custom squared! metric function
# R = :discrete (default), :continuous (leads to avoiding coincidence matrices) or precomputed (unique values of the input)
function krippendorffs_alpha(input; units::Union{Symbol,AbstractString} = :columns, metric = :nominal, R = :discrete)
    units_iterator, elementtype, insize = dispatch_units_iterator(input, Symbol(units))
    implementation, possible_responses = dispatch_responses(R, units_iterator, elementtype, insize)
    distance_metric_squared = dispatch_metric(metric, elementtype)
    return compute_alpha(implementation, possible_responses, distance_metric_squared, units_iterator)
end

const alpha = krippendorffs_alpha # enables Krippendorff.alpha
const α = krippendorffs_alpha # enable Krippendorff.α
export krippendorffs_alpha

## input type dispatch + find a meaningful elementtype + size of the whole thing

# fallback version if not specialized
# checks if Tables.jl interface implemented, 
# then if eachrow/eachcolumn work, 
# then if generic iteration works and finally fails if nothing is possible
function dispatch_units_iterator(input, units::Symbol)
    units ∈ (:rows, :columns) || error("\"units\" must be one of :rows or :columns, got: $(units)")
    userows = units === :rows # use rows as units or not (use columns then)
    if istable(input) # put this behind eachrow/eachcol? Those are probably more efficient if they exist
        # use Tables.jl interface
        return userows ? Tables.rows(input) : Tables.Columns(input)
    elseif userows && hasmethod(eachrow, Tuple{typeof(input)})
        return eachrow(input)
    elseif !userows && hasmethod(eachcol, Tuple{typeof(input)})
        return eachcol(input)
    else # assume it's generically iterable then and let it fail later if it isn't
        return input
    end
end

## dispatch possible responses

function dispatch_responses(R, itr, elementtype, insize)
    if R === :continuous 
        return (Generic(), nothing)
    elseif R === :discrete
        # get unique values of input, if there are lots of unique values 
        # in comparison to the input size, use generic version
        # also don't bother with coincidence etc if the input is tiny
        uniques = unique(Base.Iterators.flatten(itr))
        # TODO determine cutoffs
        if insize < 100 || insize / length(uniques) > 0.5
            return (Generic(), nothing)
        else
            return (Coincidence(), optimize_possible_responses(uniques, elementtype))
        end
    else # assume R is precomputed and reasonable
        return (Coincidence(), R)
    end
end

# try to find a more efficient mapping strategy than dict(unique_values => indices)
# i.e. if the input looks like a range and is sufficiently dense, then make it a range
function optimize_possible_reponses(uniquevalues, elementtype)
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
map_responses_to_indices(::Nothing) = nothing
map_responses_to_indices(::Base.OneTo) = identity
map_responses_to_indices(possible_responses::UnitRange) = Base.Fix2(-, possible_responses.start - 1)
function map_responses_to_indices(possible_responses)
    let mapping = reverse.(collect(pairs(possible_responses)))
        return x -> mapping[x]
    end
end

## dispatch inputs for metric (and assert that it fits the elementtype)
const KNOWNMETRICS = [:nominal, :interval]

function dispatch_metric(metric::Union{Symbol, AbstractString}, elementtype)
    m = Symbol(metric)
    if m === :nominal
        return (==)
    elseif m === :interval
        elementtype <: Number || @warn "Interval metric is unlikely to work with non-numeric responses."
        return (x,y)->abs2(x-y)
    else 
        raise(ArgumentError("Metric $(m) is currently not known. Try to use one of: $(KNOWNMETRICS)"))
    end
end

function dispatch_metric(metric, _) # takes everything that's not a symbol or string
    # assumes that it's a valid function defined on every pair of possible responses
    return metric
end

# from a distance metric and the set of possible responses mapping to the coincidence matrix indices, 
# compute a symmetric distance matrix once for later evaluation

function compute_distance_matrix(possible_responses, mapping_function::F, distance_metric_squared) where {F}
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
function compute_alpha(::Coincidence, possible_responses, distance_metric_squared, units_iterator)
    mapping_function = map_responses_to_indices(possible_responses)
    distances = compute_distance_matrix(possible_responses, mapping_function, distance_metric_squared)
    l = length(possible_responses)
    C = zeros(l, l) # initialize coincidence matrix
    fill_coincidence_matrix!(C, mapping_function, units_iterator)
    return compute_alpha_from_coincidence(C, distances)
end

function fill_coincidence_matrix!(C, mapping_function::F, units_iterator) where {F} # should force specialization
    for unit in units_iterator
        increment::Float64 = 1 / (length(unit) - 1)
        for i in eachindex(unit), j in eachindex(unit)
            if i<j # this assumes ordered indices, won't work for e.g. NamedTuples but I think this is rarely a problem
                cᵢ = mapping_function(unit[i])
                cⱼ = mapping_function(unit[j])
                C[ cᵢ, cⱼ ] += increment
                C[ cⱼ, cᵢ ] += increment
            elseif i===j
                cᵢ = mapping_function(unit[i])
                C[ cᵢ, cᵢ ] += 2*increment
            end
        end
    end
end

function compute_alpha_from_coincidence(C, distance_matrix)
    l = size(C,1)
    marginals = sum(C, dims=1)
    n = sum(marginals)
    disagreement_observed = sum( C[i,j] * distance_matrix[i,j] for i in 1:l for j in i+1:l)
    disagreement_expected = sum( marginals[i] * marginals[j] * distance_matrix[i,j] for i in 1:l for j in i+1:l)
    # it would be sufficient to evaluate just the upper or lower triangle of each, but I don't think this would be faster
    return 1 - (n-1)*(disagreement_observed / disagreement_expected)
end

# TODO docstring
# generic without coincidence matrix
function compute_alpha(::Generic, _, distance_metric, units_iterator)
    
end

end # module
