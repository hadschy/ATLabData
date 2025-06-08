using ATLabData
using Test

@testset "ATLabData.jl" begin
    data = load("/home/thomas/simulations/test.atlab/bPrime014200")
    println(size(data))
    data = abs(2*data - data*data + 3*data/data)^2
    data = crop(data)
    # avg = average(data)
    # avg = rms(data)
    # fig, ax, hm = visualise(data)
end
