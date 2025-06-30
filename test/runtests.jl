using ATLabData
using Test

@testset "ATLabData.jl" begin
    
    # --------------------------- ScalarData
    printstyled("Testing ScalarData loading, operations and plotting: \n", bold=true)
    @time data = load("/home/thomas/simulations/test.atlab/bPrime014200")
    @time data = abs(data*2.0f0 - data*data + 3*data/data)^2
    @time data = crop(
        data, 
        xmin = data.grid.x[round(Int, data.grid.nx/4)],
        xmax = data.grid.x[round(Int, data.grid.nx*3/4)],
        zmin = data.grid.z[round(Int, data.grid.nz/4)],
        zmax = data.grid.z[round(Int, data.grid.nz*3/4)],
    )
    @time avg = rms(data)
    @time avg = average(data)
    @time fig, ax, hm = heatmap(data)


    # --------------------------- VectorData
    printstyled("Testing VectorData loading: \n", bold=true)
    @time data = load(
        "/home/thomas/simulations/periodic/shear_layer/Re22500/RIg1.0/dU2.0/FRb0.0/VelocityVector030000.1",
        "/home/thomas/simulations/periodic/shear_layer/Re22500/RIg1.0/dU2.0/FRb0.0/VelocityVector030000.2",
        "/home/thomas/simulations/periodic/shear_layer/Re22500/RIg1.0/dU2.0/FRb0.0/VelocityVector030000.3"
    )
    # @time data = curl(data, order=2)


    # ------------------------- AveragesData
    println("Testing AveragesData loading and plotting:")
    @time data = load("/home/thomas/simulations/test.atlab/avg14200.nc", "rU")
    @time fig, ax, ln = lines(data)


    # ------------------------- From raw data
    printstyled("Testing raw file loading: \n", bold=true)
    @time data = load("/home/thomas/simulations/periodic/shear_layer/Re22500/RIg1.0/dU2.0/FRb0.0/flow.30000.1")
    @time data = abs(data*2.0f0 - data*data + 3*data/data)^2
    @time data = crop(
        data, 
        xmin = data.grid.x[round(Int, data.grid.nx/4)],
        xmax = data.grid.x[round(Int, data.grid.nx*3/4)],
        zmin = data.grid.z[round(Int, data.grid.nz/4)],
        zmax = data.grid.z[round(Int, data.grid.nz*3/4)],
    )
    @time avg = rms(data)
    @time avg = average(data)
    @time fig, ax, hm = heatmap(data)

    @time data = load(
        "/home/thomas/simulations/periodic/shear_layer/Re22500/RIg1.0/dU2.0/FRb0.0/flow.30000.1",
        "/home/thomas/simulations/periodic/shear_layer/Re22500/RIg1.0/dU2.0/FRb0.0/flow.30000.2",
        "/home/thomas/simulations/periodic/shear_layer/Re22500/RIg1.0/dU2.0/FRb0.0/flow.30000.3"
    )
end
