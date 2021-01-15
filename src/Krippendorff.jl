module Krippendorff

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
    units_iterator = dispatch_units_iterator(input, Symbol(units))
    implementation, possible_responses = dispatch_responses(R, input)
    distance_metric = dispatch_metric(metric)
    return compute_alpha(Val(implementation), possible_responses, distance_metric, units_iterator)
end

const alpha = krippendorffs_alpha # enables Krippendorff.alpha
const α = krippendorffs_alpha # enable Krippendorff.α
export krippendorffs_alpha

## input type dispatch

# fallback version if not specialized
# checks if Tables.jl interface implemented, 
# then if eachrow/eachcolumn work, 
# then if generic iteration works and finally fails if nothing is possible
function dispatch_units_iterator(input, units::Symbol)
end

## dispatch possible responses

function dispatch_responses(R, input)
    if R === :continuous 
        return (:generic, nothing)
    elseif R === :discrete
        # get unique values of input, if there are lots of unique values 
        # in comparison to the input size, use generic version
        # TODO implement that, determine cutoff
        return (:coincidence, unique(input)) # mocked for now
    else # assume R is precomputed and reasonable
        return (:coincidence, R)
    end
end

## dispatch inputs for metric



## computation backends

# TODO docstrings
# using coincidence matrix
function compute_alpha(::Val{:coincidence}, possible_responses, distance_metric, units_iterator)
end

# TODO docstring
# generic without coincidence matrix
function compute_alpha(::Val{:generic}, possible_responses, distance_metric, units_iterator)
end

end # module
