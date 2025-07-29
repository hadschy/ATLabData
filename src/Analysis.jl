module Analysis


using ..DataStructures

using Base.Threads
using Interpolations: interpolate, Gridded, Linear, Constant, NoInterp
using FiniteDifferences: central_fdm, grad


export gradient, curl


"""
    gradient(data; order=4)
Return the gradient of _data_ by using the packages _FiniteDifferences_  and 
_Iterpolations_.  
_order_ determines the numerical error order for the derivatives.
"""
# function gradient(data::ScalarData; order::Int=4)::VectorData
#     # TODO: 3D, faster?, parallise
#     println("Calculating gradient with $(nthreads()) threads...")
#     printstyled("    "*data.name*"\n", color=:cyan)
#     return gradient_multi_thread(data, order)
#     # return VectorData(
#     #     name = "gradient(" * data.name * ")",
#     #     grid = data.grid,
#     #     time = data.time,
#     #     xfield = res[1,:,:,:],
#     #     yfield = res[2,:,:,:],
#     #     zfield = res[3,:,:,:]
#     # )
# end


function gradient_single_thread!(
        res::Array{T}, data::ScalarData{T,I}, itp, order::Int, ichunk::UnitRange
    )::Array{T} where {T<:AbstractFloat, I<:Signed}
    println("Process $(threadid()) in gradient_single_thread")
    k = 1
    if ichunk[1]==1
        iminoffset = round(Int, order/2+2)
    else
        iminoffset = 0
    end
    if ichunk[end]==data.grid.nx
        imaxoffset = round(Int, order/2+2)
    else
        imaxoffset = 0
    end
    for k ∈ round(Int, order/2+2):round(Int, data.grid.nz - order/2-2)
        for j ∈ 1:data.grid.ny
            for i ∈ round(Int, ichunk[1]+iminoffset):round(Int, ichunk[end] - imaxoffset)
                buffer = grad(central_fdm(order, 1), itp, data.grid.x[i], data.grid.z[k])
                # buffer = grad(central_fdm(order, 1), itp, i, k)
                res[1,i,j,k] = buffer[1]
                res[2,i,j,k] = 0.0
                res[3,i,j,k] = buffer[2]
            end
        end
    end
    return res
end


function gradient_single_thread(
        data::ScalarData{T,I}, itp, order::Int, ichunk::UnitRange
    )::Array{T} where {T<:AbstractFloat, I<:Signed}
    println("Process $(threadid()) in gradient_single_thread")
    k = 1
    if ichunk[1]==1
        iminoffset = round(Int, order/2+2)
    else
        iminoffset = 0
    end
    if ichunk[end]==data.grid.nx
        imaxoffset = round(Int, order/2+2)
    else
        imaxoffset = 0
    end
    res = Array{T}(undef, 3, length(ichunk), data.grid.ny, data.grid.nz)
    for k ∈ round(Int, order/2+2):round(Int, data.grid.nz - order/2-2)
        for j ∈ 1:data.grid.ny
            for i ∈ round(Int, ichunk[1] + iminoffset):round(Int, ichunk[end] - imaxoffset)
                buffer = grad(central_fdm(order, 1), itp, data.grid.x[i], data.grid.z[k])
                # buffer = grad(central_fdm(order, 1), itp, i, k)
                res[1,i-ichunk[1],j,k] = buffer[1]
                res[2,i-ichunk[1],j,k] = 0.0
                res[3,i-ichunk[1],j,k] = buffer[2]
            end
        end
    end
    return res
end


function gradient(
        data::ScalarData{T,I}, order::Int
    )::VectorData{T,I} where {T<:AbstractFloat, I<:Signed}
    # TODO: 3D, faster?, parallise
    println("Calculating gradient with $(nthreads()) threads...")
    printstyled("    "*data.name*"\n", color=:cyan)
    
    res = Array{T}(undef, 3, data.grid.nx, data.grid.ny, data.grid.nz)
    k = 1
    itp = interpolate(
        (data.grid.x, data.grid.z), # Nodes of the grid
        data.field[:,k,:], # Field to be interpolated
        Gridded(Linear())
    )
    "Split the data  in chunks in the x-dimension"
    xchunks = Iterators.partition(range(1, data.grid.nx), round(Int, data.grid.nx/nthreads()))
    for xchunk in collect(xchunks)
        @sync begin
            @spawn begin
                println("I am process $(threadid()) of $(nthreads())")
                println("ichunk = $xchunk")

                # gradient_single_thread!(res, data, itp, order, xchunk)
                res[:,xchunk,:,:] .= gradient_single_thread(data, itp, order, xchunk)
                
                println("... process $(threadid()) of $(nthreads()) is done")
            end
        end
    end

    return VectorData(
        name = "gradient(" * data.name * ")",
        grid = data.grid,
        time = data.time,
        xfield = res[1,:,:,:],
        yfield = res[2,:,:,:],
        zfield = res[3,:,:,:]
    )
end


function curl(data::VectorData; order::Int=4)::VectorData
    # TODO: 3D, faster?, parallise
    
    println("Calculating curl ...")
    printstyled("    "*data.name*"\n", color=:cyan)
    
    xres = zeros(Float32, (data.grid.nx, data.grid.ny, data.grid.nz))
    yres = zeros(Float32, (data.grid.nx, data.grid.ny, data.grid.nz))
    zres = zeros(Float32, (data.grid.nx, data.grid.ny, data.grid.nz))
    
    if data.grid.nz > 1
        kmin = order÷2 + 1; kmax = data.grid.nz - order - 1
    else
        kmin = 1; kmax = data.grid.nz
    end
    if data.grid.ny > 1
        jmin = order÷2 + 1; jmax = data.grid.ny - order - 1
    else
        jmin = 1; jmax = data.grid.ny
    end
    if data.grid.nx > 1
        imin = order÷2 + 1; imax = data.grid.nx - order - 1
    else
        imin = 1; imax = data.grid.nx
    end

    """ Not the real curl: computes curl for each y-slice (2D) ... TODO """
    for j ∈ jmin:jmax
        xitp = interpolate(
            (data.grid.x, data.grid.z), # Nodes of the grid
            data.xfield[:,j,:], # Field to be interpolated
            Gridded(Linear())
        )
        # yitp = interpolate(
        #     (data.grid.x, data.grid.y), # Nodes of the grid
        #     data.yfield[:,:,k], # Field to be interpolated
        #     Gridded(Linear())
        # )
        zitp = interpolate(
            (data.grid.x, data.grid.z), # Nodes of the grid
            data.zfield[:,j,:], # Field to be interpolated
            Gridded(Linear())
        )
        for k ∈ kmin:kmax
            for i ∈ imin:imax
                xgrad = grad(central_fdm(order, 1), xitp, data.grid.x[i], data.grid.z[k])
                # ygrad = grad(central_fdm(order, 1), yitp, data.grid.x[i], data.grid.z[k])
                zgrad = grad(central_fdm(order, 1), zitp, data.grid.x[i], data.grid.z[k])

                # xres[i,j,k] = zgrad[2] - ygrad[3]
                yres[i,j,k] = xgrad[3] - zgrad[1]
                # zres[i,j,k] = ygrad[1] - xgrad[2]
            end
        end
    end
    return VectorData(
        name="curl($(data.name))", time=data.time, grid=data.grid,
        xfield=xres, yfield=yres, zfield=zres
    )
end

    
end