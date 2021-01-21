module Krippendorff

using Tables: istable, rows, Columns

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
# metric = one of (String or Symbol) nominal, ordinal, interval, squared, ratio,...? or a custom metric function
# R = :discrete (default), :continuous (leads to avoiding coincidence matrices) or precomputed (unique values of the input)
function krippendorffs_alpha(input; units::Union{Symbol,AbstractString} = :rows, metric = :nominal, R = :discrete)
    units_iterator, elementtype, insize = dispatch_units_iterator(input, Symbol(units))
    implementation, possible_responses = dispatch_responses(R, units_iterator, elementtype, insize)
    distance_metric = dispatch_metric(metric, elementtype)
    return compute_alpha(Val(implementation), possible_responses, distance_metric, units_iterator)
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
    if istable(input)
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
        return (:generic, nothing)
    elseif R === :discrete
        # get unique values of input, if there are lots of unique values 
        # in comparison to the input size, use generic version
        # also don't bother with coincidence etc if the input is tiny
        uniques = unique(Base.Iterators.flatten(itr))
        # TODO determine cutoffs
        if insize < 100 || insize / length(uniques) > 0.5
            return (:generic, nothing)
        else
            return (:coincidence, optimize_possible_responses(uniques, elementtype))
        end
    else # assume R is precomputed and reasonable
        return (:coincidence, R)
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

function dispatch_metric(metric, elementtype)

end

## computation backends

# TODO docstrings
# using coincidence matrix
function compute_alpha(::Val{:coincidence}, possible_responses, distance_metric, units_iterator)
    mapping_function = map_responses_to_indices(possible_responses)
    l = length(possible_responses)
    C = zeros(Int, (l, l)) # initialize coincidence matrix
    fill_coincidence_matrix!(C, mapping_function, units_iterator)
    return compute_alpha_from_coincidence(C, distance_metric)
end

function fill_coincidence_matrix!(C, mapping_function, units_iterator)
    for unit in units_iterator
#TODO
        for i in eachindex(unit), j in eachindex(unit)
            if i<j
                
            elseif i===j

            end
        end

    end
end

function compute_alpha_from_coincidence(C, distance_metric)
end

# TODO docstring
# generic without coincidence matrix
function compute_alpha(::Val{:generic}, _, distance_metric, units_iterator)
end

end # module
