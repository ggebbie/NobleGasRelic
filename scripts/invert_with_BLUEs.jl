    @dim Vintage "vintage"
    @dim InteriorLocation "interior location"

    yr = u"yr"
    nτ = 5 # how much of a lag is possible?
    lags = (0:(nτ-1))yr
    surfaceregions = [:NATL,:ANT,:SUBANT]
    years = (1990:2000)yr
    n = length(surfaceregions)

    M = source_water_matrix_with_lag(surfaceregions,lags,years)

    x= source_water_solution(surfaceregions,years)

    # Run model to predict interior location temperature
    # convolve E and x
    y = convolve(x,M)

    # could also use this format
    y = predictobs(convolve,x,M)

    ## invert for y for x̃
    # Given, M and y. Make first guess for x.x
    # add adjustment
    # DimArray is good enough. This is an array, not necessarily a matrix.
    x₀ = DimArray(zeros(size(x))K,(Ti(years),last(dims(M))))

    # probe to get E matrix. Use function convolve
    E = impulseresponse(convolve,x₀,M)
    
    @test (E*UnitfulMatrix(vec(x₀)))[1] .== ustrip(convolve(x₀, M))
    # Does E matrix work properly?
    ỹ = E*UnitfulMatrix(vec(x))
    @test isapprox(y,getindexqty(ỹ,1))

    x̂ = E\UnitfulMatrix([y]) 
    @test isapprox(y,getindexqty(E*x̂,1))

    # now in a position to use BLUEs to solve
    # should handle matrix left divide with
    # unitful scalars in UnitfulLinearAlgebra
    
    σₙ = 0.01
    σₓ = 100.0

    Cnn = UnitfulMatrix(Diagonal(fill(σₙ^2,length(y))),fill(unit.(y).^1,length(y)),fill(unit.(y).^-1,length(y)),exact=false)

    Cxx = UnitfulMatrix(Diagonal(fill(σₓ^2,length(x₀))),unit.(x₀)[:],unit.(x₀)[:].^-1,exact=false)

    problem = UnderdeterminedProblem(UnitfulMatrix([y]),E,Cnn,Cxx,x₀)
    x̃ = solve(problem)
    @test within(y[1],getindexqty(E*x̃.v,1),3σₙ) # within 3-sigma

    @test cost(x̃,problem) < 5e-2 # no noise in ob
end
