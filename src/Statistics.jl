module Statistics


import StatsBase: mean

using ..DataStructures

export average, rms, mean, flucs


"""
    average(data) -> AveragesData
Compute the arithmetric mean along the second dimension while preserving metadata.

    average(data; coord::Int=...) -> AveragesData
Compute the arithmetric mean along the dimension given by `coord` as Int.
"""
function average(data::ScalarData; coord=3)::AveragesData
    println("Computing averages along coord=$coord ...")
    printstyled("   $(data.name)  \n")
    if coord==1
        res = zeros(eltype(data)[1], data.grid.nx)
        for i ∈ 1:data.grid.nx
            res[i] = sum(view(data.field, i, :, :)) / (data.grid.ny*data.grid.nz)
        end
        range = data.grid.x
    elseif coord==2
        res = zeros(eltype(data)[1], data.grid.ny)
        for j ∈ 1:data.grid.ny
            res[j] = sum(view(data.field, :, j, :)) / (data.grid.nx*data.grid.nz)
        end
        range = data.grid.y
    elseif coord==3
        res = zeros(eltype(data)[1], data.grid.nz)
        for k ∈ 1:data.grid.nz
            res[k] = sum(view(data.field, :, :, k)) / (data.grid.nx*data.grid.ny)
        end
        range = data.grid.z
    else
        error("coord has be in {1,2,3}")
    end
    name = "avg$coord($(data.name))"
    return AveragesData(name=name, time=data.time, range=range, field=res)
end


function rms(data::ScalarData; coord=3)::AveragesData
    println("Computing rms along coord=$coord ...")
    printstyled("   $(data.name) \n")
    # nv = size(data)[coord]
    # nh = data.grid.nx*data.grid.ny*data.grid.nz / nv
    if coord==1
        res = zeros(eltype(data)[1], data.grid.nx)
        for i ∈ 1:data.grid.nx
            res[i] = sum(view(data.field, i, :, :).^2)
            res[i] = sqrt(res[i]/(data.grid.ny*data.grid.nz))
        end
        res = data.grid.x
    elseif coord==2
        res = zeros(eltype(data)[1], data.grid.ny)
        for i ∈ 1:data.grid.ny
            res[i] = sum(view(data.field, :, i, :).^2)
            res[i] = sqrt(res[i]/(data.grid.nx*data.grid.nz))
        end
        range = data.grid.y
    elseif coord==3
        res = zeros(eltype(data)[1], data.grid.nz)
        for i ∈ 1:data.grid.nz
            res[i] = sum(view(data.field, :, :, i).^2)
            res[i] = sqrt(res[i]/(data.grid.nx*data.grid.nz))
        end
        range = data.grid.z
     else
        error("Give coord as Int in {1,2,3}.")
    end
    name = "rms$(coord)($(data.name))"
    return AveragesData(name=name, time=data.time, range=range, field=res)
end


"""
Computes the mean field and returns it in same dimensions as data.
"""
function mean(data::VectorData)::VectorData
    resx = zeros(eltype(data.xfield), size(data.xfield))
    resy = zeros(eltype(data.yfield), size(data.yfield))
    resz = zeros(eltype(data.zfield), size(data.zfield))
    for k ∈ 1:data.grid.nz
        resx[:,:,k] .= sum(data.xfield[:,:,k])./(data.grid.nx*data.grid.ny)
        resy[:,:,k] .= sum(data.yfield[:,:,k])./(data.grid.nx*data.grid.ny)
        resz[:,:,k] .= sum(data.zfield[:,:,k])./(data.grid.nx*data.grid.ny)
    end
    return VectorData("mean($(data.name))", data.grid, data.time, resx, resy, resz)
end


function mean(data::ScalarData)::ScalarData
    res = zeros(eltype(data.field), size(data.field))
    for k ∈ 1:data.grid.nz
        res[:,:,k] .= sum(data.field[:,:,k])./(data.grid.nx*data.grid.ny)
    end
    return ScalarData("mean($(data.name))", data.grid, data.time, res)
end


flucs(data::VectorData)::VectorData = data - mean(data)


end