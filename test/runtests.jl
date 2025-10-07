using ATLabData
using Test

@testset "ATLabData.jl" begin
    
    grid = loadgrid("/home/thomas/simulations/tmp/WTD2025/dU2.0/grid")
    data = init(grid)
    load!(data, "/home/thomas/simulations/tmp/WTD2025/dU2.0/VorticityVector019200.2")
end
