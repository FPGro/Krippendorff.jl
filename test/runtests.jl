using Test, Krippendorff, Tables

# write some tests
# should do some unit tests and then some tests against known reliability values

# krippendorffs_alpha

# examples from: "Computing Krippendorff’s Alpha-Reliability, Klaus Krippendorff, kkrippendorff at asc.upenn.edu, 2011.1.25"
binarytest = [ 0 1 0 0 0 0 0 0 1 0   # α ≈ 0.095 in his example
               1 1 1 0 0 1 0 0 0 0 ] 
@test round(Krippendorff.alpha(binarytest, units = :cols), sigdigits = 2) == 0.095 
@test krippendorffs_alpha(binarytest', units = :rows) == Krippendorff.alpha(binarytest, units = :columns)
# the next test passes only with approx, indicating some minor difference between the approaches
@test Krippendorff.alpha(binarytest, R = 0:1) ≈ Krippendorff.alpha(binarytest, R = :continuous)

# _dispatch_units_iterator
# _units_eltype_and_size
# _units_eltype_and_size(Matrix)
# _dispatch_responses(precomputed R)
# _dispatch_responses(Symbol/String)
# _optimize_possible_reponses
# _map_responses_to_indices(mapping)
# _dispatch_metric(Symbol/String)
@test Krippendorff._dispatch_metric("nominal", nothing) == (!=)
kripp_intervall = Krippendorff._dispatch_metric(:interval, Int) 
test_intervall = (x,y)->abs2(x-y)
intervaltest_rands = rand(Int,(100,2))
@test all([kripp_intervall(rands[1],rands[2]) for rands in eachrow(intervaltest_rands)] .== [test_intervall(rands[1],rands[2]) for rands in eachrow(intervaltest_rands)])
@test_throws Krippendorff._dispatch_metric(:blob, Int) ErrorException
testmetric = (x,y)->(x-y)
@test Krippendorff._dispatch_metric(testmetric, Int) == testmetric
# _compute_distance_matrix
# _compute_alpha
# _fill_coincidence_matrix
# _evaluate_coincidence_matrix
# compute_alpha_generical
# compute_alpha_with_coincidences
# prepare_iterator
@test all(iterate(Krippendorff.prepare_iterator(reshape(1:20,(4,5))))[1] .== [1,5,9,13,17])
@test collect(Iterators.flatten(Krippendorff.prepare_iterator([[1,2,missing],[missing,missing,missing],[3,missing]]))) == [1,2,3]
# istable