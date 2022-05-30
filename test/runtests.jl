using Test
using Krippendorff
using Tables

## krippendorffs_alpha

# examples from: "Computing Krippendorff’s Alpha-Reliability, Klaus Krippendorff, kkrippendorff at asc.upenn.edu, 2011.1.25"
let binarytest = [ 0 1 0 0 0 0 0 0 1 0   # α ≈ 0.095 in his example
                   1 1 1 0 0 1 0 0 0 0 ] 
    @test round(Krippendorff.alpha(binarytest, units = :cols), sigdigits = 2) == 0.095 
    @test krippendorffs_alpha(binarytest', units = :rows) == Krippendorff.alpha(binarytest, units = :columns)
# the next test passes only with approx, indicating some minor difference between the approaches
    @test Krippendorff.alpha(binarytest, R = 0:1) ≈ Krippendorff.alpha(binarytest, R = :continuous)
end

## _dispatch_units_iterator
@test_throws ArgumentError Krippendorff._dispatch_units_iterator([1,2,3], :foo)
@test_throws ArgumentError Krippendorff._dispatch_units_iterator([], :rows) # empty input
let input1D = [[1,2,3],[4,5,5],[2,2,2]]
    rowitr, rowtype, rowsize = Krippendorff._dispatch_units_iterator(input1D, :rows)
    colitr, coltype, colsize = Krippendorff._dispatch_units_iterator(input1D, :cols)
    @test rowitr == colitr
    @test rowtype == coltype == eltype(eltype(input1D))
    @test rowsize == colsize == length(collect(Iterators.flatten(input1D)))
end
let input2D = [[1,2,3] [4,5,5] [2,2,2]]
    rowitr, rowtype, rowsize = Krippendorff._dispatch_units_iterator(input2D, :rows)
    colitr, coltype, colsize = Krippendorff._dispatch_units_iterator(input2D, :cols)
    @test rowitr != colitr
    @test rowtype == coltype == eltype(input2D)
    @test rowsize == colsize == length(input2D)
end

## _units_eltype_and_size

## _units_eltype_and_size(Matrix)

## _dispatch_responses(precomputed R)

## _dispatch_responses(Symbol/String)

## _optimize_possible_reponses

## _map_responses_to_indices(mapping)

## _dispatch_metric(Symbol/String)
@test Krippendorff._dispatch_metric("nominal", nothing) == (!=)
let kripp_intervall = Krippendorff._dispatch_metric(:interval, Int),
    test_intervall = (x,y)->abs2(x-y),
    intervaltest_rands = rand(Int,(100,2))
    
    @test all( [kripp_intervall(rands[1],rands[2]) for rands in eachrow(intervaltest_rands)] .== 
                [test_intervall(rands[1],rands[2]) for rands in eachrow(intervaltest_rands)])
end
@test_throws ArgumentError Krippendorff._dispatch_metric(:blob, Int)
let testmetric = (x,y)->(x-y)
    @test Krippendorff._dispatch_metric(testmetric, Int) == testmetric
end

## _compute_distance_matrix

## _compute_alpha

## _fill_coincidence_matrix

## _evaluate_coincidence_matrix

## compute_alpha_generical

## compute_alpha_with_coincidences

## prepare_iterator
@test all(iterate(Krippendorff.prepare_iterator(reshape(1:20,(4,5))))[1] .== [1,5,9,13,17])
@test collect(Iterators.flatten(Krippendorff.prepare_iterator([[1,2,missing],[missing,missing,missing],[3,missing]]))) == [1,2,3]

## istable
# not strictly necessary, because this only wraps Tables.istable + fancy output
let mat = rand((10,20))
    @test Krippendorff.istable(mat) === false
    @test Krippendorff.istable(Tables.table(mat)) === true
end
for reallyatable in [   (a = [1,2,3], b = [4.,5.,6.]) # NamedTuple of Vectors
                        Dict([:a => [1,2,3], :b => [4.,5.,6.]]) # Dict{Symbol,Vector}
                        [(a = 1,b = 2), (a = 3, c = 4)] # Vector of NamedTuples
                        [Dict([:a => 2, :b => 4]), Dict([:a => 01, :c => 10])] # Vector of Dict{Symbol,_}
                    ]
    @test Krippendorff.istable(reallyatable) === true
end
for notatable in [  (a = [1,2,3], b = nothing) # NamedTuple of Any
                    Dict(["a" => [1,2,3], "b" => [4.,5.,6.]]) # Dict{not_Symbol,Vector}
                    Dict([:a => [1,2,3], :b => nothing]) # Dict{Symbol,not_Vector}
                    [Dict([1 => 2, 3 => 4]), Dict([1 => 01, 2 => 10])] # Vector of Dict{not_Symbol,_}
                ]
    @test Krippendorff.istable(notatable) === false
end
