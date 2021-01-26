module Krippendorff

using Tables

# TODO:
# - documentation
# - check for correctness
# - add remaining metrices
# - add tests

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
    units_iterator, elementtype, inputsize = dispatch_units_iterator(input, Symbol(units))
    implementation, possible_responses = dispatch_responses(R, units_iterator, elementtype, inputsize)
    distance_metric_squared = dispatch_metric(metric, elementtype)
    return compute_alpha(implementation, possible_responses, distance_metric_squared, units_iterator)
end

const alpha = krippendorffs_alpha # enables Krippendorff.alpha
export krippendorffs_alpha

## input type dispatch + find a meaningful elementtype + size of the whole thing

# fallback version if not specialized
# checks if Tables.jl interface implemented, 
# then if eachrow/eachcolumn work, 
# then if generic iteration works and finally fails if nothing is possible
function dispatch_units_iterator(input, units::Symbol)
    units ∈ (:rows, :columns, :cols) || error("\"units\" must be one of :rows or :col(umn)s, got: $(units)")
    userows = units === :rows # use rows as units or not (use columns then)
    units, elementtype, sizeestimate = units_eltype_and_size(input, userows)
    # if missings are allowed, wrap every unit in a skipmissings
    prepared_iter, prepared_eltype = Missing<:elementtype ? ((skipmissing(u) for u in units), nonmissingtype(elementtype)) : (units, elementtype)
    return (prepared_iter, prepared_eltype, sizeestimate)
end

# use multiple dispatch to find an appropriate units_iterator, eltype and size for different types

function units_eltype_and_size(input, userows::Bool)
    # fallback, check Tables.jl interface via istable, else assume generic iterability 
    if Tables.istable(input)
        if userows
            rows = Tables.rows(input)
            schema = Tables.schema(rows)
            if isnothing(schema)
                # infer from first row
                @info "Schema-less table detected, infering layout from first row only. This may lead to unexpected results."
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
                @info "Schema-less table detected, infering types from first column only. This may lead to unexpected results."
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
        @info "Input has no notion of rows and columns, iterating it is expected to yield units."
        return (units, elementtype, sizeestimate)
    end
end

function units_eltype_and_size(input::AbstractMatrix, userows::Bool)
    iter = userows ? eachrow(input) : eachcol(input)
    return (iter, eltype(input), *(length(axes(input,1)), length(axes(input,2))))
end

## dispatch possible responses

function dispatch_responses(R, _,_,_)
    # if not a Symbol: assume R is precomputed and reasonable
    return (Coincidence(), R)
end

function dispatch_responses(R::Union{Symbol,AbstractString}, itr, elementtype, inputsize)
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
    let mapping = Dict(reverse.(collect(pairs(possible_responses))))
        return x -> mapping[x]
    end
end

## dispatch inputs for metric (and assert that it fits the elementtype)
const KNOWNMETRICS = [:nominal, :interval]

function dispatch_metric(metric::Union{Symbol, AbstractString}, elementtype)
    m = Symbol(metric)
    if m === :nominal
        return (!=)
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

function compute_alpha_from_coincidence(C, distance_matrix)
    l = size(C,1)
    marginals = sum(C, dims=1)
    n = sum(marginals)
    disagreement_observed = sum( C[i,j] * distance_matrix[i,j] for i in 1:l for j in i+1:l)
    disagreement_expected = sum( marginals[i] * marginals[j] * distance_matrix[i,j] for i in 1:l for j in i+1:l)
    # it would be sufficient to evaluate just the upper or lower triangle of each, but I don't think this would be faster
    return 1. - (n-1)*(disagreement_observed / disagreement_expected)
end

# TODO docstring
# generic without coincidence matrix
function compute_alpha(::Generic, _, distance_metric, units_iterator)
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

end # module
