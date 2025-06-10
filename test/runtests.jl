using ATLabData
using Test

@testset "ATLabData.jl" begin
    @time data = load("/home/thomas/simulations/test.atlab/bPrime014200")
    println(size(data))
    @time data = abs(data*2.0f0 - data*data + 3*data/data)^2
    @time data = crop(
        data, 
        xmin = data.grid.x[round(Int, data.grid.nx/4)],
        xmax = data.grid.x[round(Int, data.grid.nx*3/4)],
        zmin = data.grid.z[round(Int, data.grid.nz/4)],
        zmax = data.grid.z[round(Int, data.grid.nz*3/4)],
    )
    
    @time avg = average(data)
    @time avg = rms(data)
    @time fig, ax, hm = heatmap(data)
    @time data = load("/home/thomas/simulations/test.atlab/avg14200.nc", "rU")
    @time fig, ax, ln = lines(data)
end
