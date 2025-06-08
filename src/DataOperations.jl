module DataOperations

import Base: +, -, *, /, ^, size, display, abs, convert, deepcopy, eltype

using Base.Threads
using Interpolations: interpolate, Gridded, Linear
using FiniteDifferences: central_fdm, grad

using ..DataStructures

export size, display, +, -, *, /, ^, abs, log, view, convert, eltype
export crop, gradient, norm, curl, logarithm
export component
export rms, average


# ------------------------------------------------------------------------------
# --------------- Base function methods for composite types --------------------
# ------------------------------------------------------------------------------

"""
    size(data)
Returns `sìze(data.field)`.
"""
size(data::ScalarData)::Tuple = size(data.field)
size(data::VectorData)::Tuple = (size(data.xfield), size(data.yfield), size(data.zfield))
size(data::AveragesData)::Tuple = size(data.field)


"""
    data1 - data2
Subtract the _field_ of _data2_ from the _field_ of _data1_. 
Does the same as '''subtract(data1, data2)'''
"""
-(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    name = "$(data1.name)-$(data2.name)",
    field = data1.field .- data2.field,
    time = data1.time,
    grid = data1.grid
)
-(data1::VectorData, data2::VectorData)::VectorData = VectorData(
    name = "$(data1.name)-$(data2.name)",
    grid = data1.grid,
    time = data1.time,
    xfield = data1.xfield .- data2.xfield,
    yfield = data1.yfield .- data2.yfield,
    zfield = data1.zfield .- data2.zfield
)
-(data1::AveragesData, data2::AveragesData)::AveragesData = AveragesData(
    name = "$(data1.name)-$(data2.name)",
    time = data1.time,
    range = data1.range,
    field = data1.field .- data2.field
)


"""
    data1 + data2
Add the _fields_ of _data1_ and _data2_. Does the same as '''add(data1, data2)'''
"""
+(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    field = data1.field .+ data2.field,
    time = data1.time,
    grid = data1.grid,
    name = data1.name*"+"*data2.name
)
+(data1::VectorData, data2::VectorData)::VectorData = VectorData(
    name = data1.name * "+" * data2.name,
    grid = data1.grid,
    time = data1.time,
    xfield = data1.xfield .+ data2.xfield,
    yfield = data1.yfield .+ data2.yfield,
    zfield = data1.zfield .+ data2.zfield
)
+(data1::AveragesData, data2::AveragesData)::AveragesData = AveragesData(
    name = "$(data1.name)+$(data2.name)",
    time = data1.time,
    range = data1.range,
    field = data1.field .+ data2.field
)


"""
    data * factor
Multiply the _field_ of _data_ with _factor_. Does the same as _rescale_.

    data1 * data2
Vectorized multiplication of the fields of _data1_ and _data2_
"""
*(data::ScalarData, factor::Real)::ScalarData = ScalarData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    field = data.field .* convert(Float32, factor)
)
*(data::VectorData, factor::Real)::VectorData = VectorData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    xfield = data.xfield * convert(Float32, factor),
    yfield = data.yfield * convert(Float32, factor),
    zfield = data.zfield * convert(Float32, factor)
)
*(data::AveragesData, factor::Real)::AveragesData = AveragesData(
    name = data.name,
    time = data.time,
    range = data.range,
    field = data.field .* convert(Float32, factor)
)
*(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    grid = data1.grid,
    field = data1.field .* data2.field
)
*(data1::AveragesData, data2::AveragesData) = AveragesData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    range = data1.range,
    field = data1.field .* data2.field
)
*(factor::Real, data::AbstractData)::AbstractData = data*factor


"""
    data / number
Divide the _field_ of _data_ with _number__.

    data1 / data2
Vectorized division of the fields of _data1_ and _data2_.
"""
/(data::ScalarData, factor::Real)::ScalarData = ScalarData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    field = data.field ./ convert(Float32, factor)
)
/(data::VectorData, factor::Real)::VectorData = VectorData(
    name = data.name,
    time = data.time,
    grid = data.grid,
    xfield = data.xfield ./ convert(Float32, factor),
    yfield = data.yfield ./ convert(Float32, factor),
    zfield = data.zfield ./ convert(Float32, factor)
)
/(data::AveragesData, factor::Real)::AveragesData = AveragesData(
    name = data.name,
    time = data.time,
    range = data.range,
    field = data.field ./ convert(Float32, factor)
)
/(data1::ScalarData, data2::ScalarData)::ScalarData = ScalarData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    grid = data1.grid,
    field = data1.field ./ data2.field
)
/(data1::AveragesData, data2::AveragesData) = AveragesData(
    name = "$(data1.name)*$(data2.name)",
    time = data1.time,
    range = data1.grid,
    field = data1.field ./ data2.field
)


"""
    data^exponent
Exponentiation of the _Data_ type. Returns _Data_  with data.field.^2 while 
maintaining the metadata.
"""
^(data::ScalarData, exponent::Real)::ScalarData = ScalarData(
    name = data.name, 
    time = data.time, 
    grid = data.grid,
    field = data.field.^exponent
)
^(data::AveragesData, exponent::Real)::AveragesData = AveragesData(
    name = data.name, 
    time = data.time, 
    range = data.range,
    field = data.field.^exponent
)


"""
    abs(data)
Returns absolute values of _data_ while preserving the metadata. 
For _VectorData_ the euclidian norm is return.
"""
abs(data::ScalarData) = ScalarData(
    field = abs.(data.field),
    grid = data.grid,
    name = "abs(" * data.name * ")",
    time = data.time
)
abs(data::VectorData) = norm(data)
norm(data::VectorData)::ScalarData = ScalarData(
    field = sqrt.(data.xfield.^2 .+ data.yfield.^2 .+ data.zfield.^2),
    grid = data.grid,
    name = "abs(" * data.name * ")",
    time = data.time
)


"""
    log(data)
Returns the logarithm while preserving the metadata. Where the result is Inf 
    or -Inf the value is replaced by zero.
"""
log(b::Base.Complex, data::ScalarData) = ScalarData(
    name = "log("*data.name*")",
    time = data.time,
    grid = data.grid,
    field = replace(log.(b, data.field), Inf=>0.0, -Inf=>0.0)
)
log(data::ScalarData) = ScalarData(
    name = "log("*data.name*")",
    time = data.time,
    grid = data.grid,
    field = replace(log.(data.field), Inf=>0.0, -Inf=>0.0)
)


"""
    convert(T, data)
Convert the entries in data to T.
"""
convert(T::Type{<:AbstractFloat}, grid::Grid)::Grid{T,Int32} = Grid{T,Int32}(
    grid.nx, grid.ny, grid.nz, 
    grid.scalex, grid.scaley, grid.scalez, 
    grid.x, grid.y, grid.z
)
convert(T::Type{<:AbstractFloat}, data::ScalarData)::ScalarData{T,Int32} = ScalarData{T,Int32}(
    data.name, convert(T, data.grid), data.time, data.field
)


"""
    eltype(data)
"""
eltype(grid::Grid)::Tuple = (eltype(grid.x), eltype(grid.nx))
eltype(data::ScalarData)::Tuple = eltype(data.grid)


# """
#     deepcopy(data, newdata)
# Copy _data_ into _newdata_.
# """
# # deepcopy(a::ScalarData, b::ScalarData) = 



"""
    display(data)
Show the available atributes of _data_.
"""
function display(data::Grid)
    println(typeof(data), " with attributes: ")
    print("   nx: "); println(data.nx)
    print("   ny: "); println(data.ny)
    print("   nz: "); println(data.nz)
    print("   scalex: "); println(data.scalex)
    print("   scaley: "); println(data.scaley)
    print("   scalez: "); println(data.scalez)
    print("   x: "); println(typeof(data.x))
    print("   y: "); println(typeof(data.y))
    print("   z: "); println(typeof(data.z))
end
function display(data::ScalarData)
    println(typeof(data), " with attributes: ")
    print("   name: "); println(data.name)
    print("   grid: "); println(typeof(data.grid))
    print("   field: "); println(typeof(data.field))
    print("   time: "); println(data.time)
    println("   Use display(data.grid) and display(data.field) to print the corresponding attributes")
end
function display(data::VectorData)
    println(typeof(data), " with attributes: ")
    print("   name: "); println(data.name)
    print("   grid: "); println(typeof(data.grid))
    print("   xfield: "); println(typeof(data.xfield))
    print("   yfield: "); println(typeof(data.yfield))
    print("   zfield: "); println(typeof(data.zfield))
    print("   time: "); println(data.time)
    println("   Use display(data.grid) and display(data.*field) to print the corresponding attributes")
end
function display(data::AveragesData)
    println(typeof(data), " with attributes: ")
    print("   name: "); println(data.name)
    print("   range: "); println(typeof(data.range))
    print("   field: "); println(typeof(data.field))
    print("   time: "); println(data.time)
end


# ------------------------------------------------------------------------------
# --------------------- Differential operations --------------------------------
# ------------------------------------------------------------------------------
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


function gradient_single_thread!(res::Array{Float32}, data::ScalarData, itp, order::Int, ichunk::UnitRange)::Array{Float32}
    println("Process $(threadid()) in gradient_single_thread")
    k = 1
    if ichunk[1]==1
        iminoffset = round(Int, order/2+1)
    else
        iminoffset = 0
    end
    if ichunk[end]==data.grid.nx
        imaxoffset = round(Int, order/2+1)
    else
        imaxoffset = 0
    end
    for j ∈ round(Int, order/2+1):round(Int, data.grid.ny - order/2-1)
        for i ∈ round(Int, ichunk[1]+iminoffset):round(Int, ichunk[end] - imaxoffset)
            buffer = grad(central_fdm(order, 1), itp, data.grid.x[i], data.grid.y[j])
            res[1,i,j,k] = buffer[1]
            res[2,i,j,k] = buffer[2]
            res[3,i,j,k] = 0.0
        end
    end
    return res
end


function gradient(data::ScalarData, order::Int)::VectorData
    # TODO: 3D, faster?, parallise
    println("Calculating gradient with $(nthreads()) threads...")
    printstyled("    "*data.name*"\n", color=:cyan)
    
    res = Array{Float32}(undef, 3, data.grid.nx, data.grid.ny, data.grid.nz)
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

                gradient_single_thread!(res, data, itp, order, xchunk)
                
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
        kmin = round(Int, order/2+1); kmax = round(Int, data.grid.nz-order/2-1)
    else
        kmin = 1; kmax = data.grid.nz
    end
    if data.grid.ny > 1
        jmin = round(Int, order/2+1); jmax = round(Int, data.grid.ny-order/2-1)
    else
        jmin = 1; jmax = data.grid.ny
    end
    if data.grid.nx > 1
        imin = round(Int, order/2+1); imax = round(Int, data.grid.nx-order/2-1)
    else
        imin = 1; imax = data.grid.nx
    end

    """ Not the real curl: computes curl for each z-slice ... TODO """
    for k ∈ kmin:kmax
        xitp = interpolate(
            (data.grid.x, data.grid.y), # Nodes of the grid
            data.xfield[:,:,k], # Field to be interpolated
            Gridded(Linear())
        )
        yitp = interpolate(
            (data.grid.x, data.grid.y), # Nodes of the grid
            data.yfield[:,:,k], # Field to be interpolated
            Gridded(Linear())
        )
        # zitp = interpolate(
        #     (data.grid.x, data.grid.y), # Nodes of the grid
        #     data.zfield[:,:,k], # Field to be interpolated
        #     Gridded(Linear())
        # )
        for j ∈ jmin:jmax
            for i ∈ imin:imax
                xgrad = grad(central_fdm(order, 1), xitp, data.grid.x[i], data.grid.y[j])
                ygrad = grad(central_fdm(order, 1), yitp, data.grid.x[i], data.grid.y[j])
                # zgrad = grad(central_fdm(order, 1), zitp, data.grid.x[i], data.grid.y[j])

                # xres[i,j,k] = zgrad[2] - ygrad[3]
                # yres[i,j,k] = xgrad[3] - zgrad[1]
                zres[i,j,k] = ygrad[1] - xgrad[2]
            end
        end
    end
    return VectorData(
        name="curl($(data.name))", time=data.time, grid=data.grid,
        xfield=xres, yfield=yres, zfield=zres
    )
end


# ------------------------------------------------------------------------------
# ------------------- Statistical operations -----------------------------------
# ------------------------------------------------------------------------------
"""
    average(data) -> AveragesData
Compute the arithmetric mean along the second dimension while preserving metadata.

    average(data; coord::Int=...) -> AveragesData
Compute the arithmetric mean along the dimension given by `coord` as Int.
"""
function average(data::ScalarData; coord=3)::AveragesData
    if coord==1
        res = zeros(Float32, data.grid.nx)
        for i ∈ 1:data.grid.nx
            res[i] = sum(data.field[i,:,:]) / (data.grid.ny*data.grid.nz)
        end
        range = data.grid.x
        name = "avg1(" * data.name * ")"
    elseif coord==2
        res = zeros(Float32, data.grid.ny)
        for j ∈ 1:data.grid.ny
            res[j] = sum(data.field[:,j,:]) / (data.grid.nx*data.grid.nz)
        end
        range = convert(Vector{Float32}, data.grid.y)
        name = "avg2(" * data.name * ")"
    elseif coord==3
        res = zeros(Float32, data.grid.nz)
        for k ∈ 1:data.grid.nz
            res[k] = sum(data.field[:,:,k]) / (data.grid.nx*data.grid.ny)
        end
        range = convert(Vector{Float32}, data.grid.z)
        name = "avg3(" * data.name * ")"
    else
        error("coord has be in {1,2,3}")
    end
    return AveragesData(name=name, time=data.time, range=range, field=res)
end


function rms(data::ScalarData; coord=3)::AveragesData
    nv = size(data)[coord]
    nh = data.grid.nx*data.grid.ny*data.grid.nz / nv
    res = zeros(nv)
    if coord==1
        for i ∈ 1:nv
            s = sum(view(data.field, i, :, :).^2)
            res[i] = sqrt(s/nh)
        end
        range = convert(Vector{Float32}, data.grid.x)
    elseif coord==2
        for i ∈ 1:nv
            s = sum(view(data.field, :, i, :).^2)
            res[i] = sqrt(s/nh)
        end
        range = convert(Vector{Float32}, data.grid.y)
    elseif coord==3
        for i ∈ 1:nv
            s = sum(view(data.field, :, :, i).^2)
            res[i] = sqrt(s/nh)
        end
        range = convert(Vector{Float32}, data.grid.z)
     else
        error("Give coord as Int in {1,2,3}.")
    end
    return AveragesData(
        name = "rms$(coord)($(data.name))", 
        time = data.time,
        range = range,
        field = convert(Vector{Float32}, res)
    )
end


# ------------------------------------------------------------------------------
# -------------------- Additional operations -----------------------------------
# ------------------------------------------------------------------------------
"""
    component(data, field) -> ScalarData
Extract from data of type _VectorData_ the component _field_ 
while keeping the metadata.
"""
function component(data::VectorData, field::Symbol)::ScalarData 
    if field == :x
        return ScalarData(
            name = data.name * " - "*string(field),
            time = data.time,
            grid = data.grid,
            field = data.xfield
        )
    elseif field == :y
        return ScalarData(
            name = data.name * " - "*string(field),
            time = data.time,
            grid = data.grid,
            field = data.yfield
        )
    elseif field == :z
        return ScalarData(
            name = data.name * " - "*string(field),
            time = data.time,
            grid = data.grid,
            field = data.zfield
        )
    else
        error("Choose the x, y, or z component of the VectorData")
    end
end


function crop(
        data::ScalarData; 
        xmin = data.grid.x[1], xmax = data.grid.x[end], 
        ymin = data.grid.y[1], ymax = data.grid.y[end],
        zmin = data.grid.z[1], zmax = data.grid.z[end]
    )::ScalarData
    println("Croping ...")
    printstyled("   "*data.name, "\n", color=:cyan)
    kmin::Int=1; kmax::Int=data.grid.nz
    jmin::Int=1; jmax::Int=data.grid.ny
    imin::Int=1; imax::Int=data.grid.nx
    inrange = false; was_inrange = false
    for k in eachindex(data.field[1,1,:])
        z = data.grid.z[k]
        if zmin <= z <= zmax
            inrange = true
        else
            inrange = false
        end
        if inrange && !(was_inrange)
            kmin = k
        end
        if !(inrange) && was_inrange
            kmax = k-1
        end
        if inrange
            was_inrange = true
        else
            was_inrange = false
        end
    end
    inrange = false; was_inrange = false
    for j in eachindex(data.field[1,:,1])
        y = data.grid.y[j]
        if ymin <= y <= ymax
            inrange = true
        else
            inrange = false
        end
        if inrange && !(was_inrange)
            jmin = j
        end
        if !(inrange) && was_inrange
            jmax = j-1
        end
        if inrange
            was_inrange = true
        else
            was_inrange = false
        end
    end
    inrange = false; was_inrange = false
    for i in eachindex(data.field[:,1,1])
        x = data.grid.x[i]
        if xmin <= x <= xmax
            inrange = true
        else
            inrange = false
        end
        if inrange && !(was_inrange)
            imin = i
        end
        if !(inrange) && was_inrange
            imax = i-1
        end
        if inrange
            was_inrange = true
        else
            was_inrange = false
        end
    end
    return ScalarData{eltype(data)[1], eltype(data)[2]}(
        data.name*"-crop",
        Grid{eltype(data)[1], eltype(data)[2]}(
            imax + 1 - imin,
            jmax + 1 - jmin,
            kmax + 1 - kmin,
            xmax - xmin,
            ymax - ymin,
            zmax - zmin,
            data.grid.x[imin:imax],
            data.grid.y[jmin:jmax],
            data.grid.z[kmin:kmax]
        ), 
        data.time,
        data.field[imin:imax,jmin:jmax,kmin:kmax],
    )
end


#-------------------------------------------------------------------------------
#                           Under Construction
#-------------------------------------------------------------------------------
# TODO See DataStructures.jl

# """
#     view(data, i, j, k)
# Method for the _Base_ function ```view``` for the _Data_ type. 
# _i, j, k_ can be either particular indices or ranges.

# _Examples_:
# ```
#     view(data, 2, 934, 1)
#     view(data, 1:50, 45:3, 4:5)
# ```
# """
# view(data::ScalarData, i::Int, j::Int, k::Int) = ScalarData(
#     name = data.name,
#     time = data.time,
#     grid = data.grid,
#     field = view(data.field, i, j, k)
# )
# view(data::ScalarData, irange::UnitRange, jrange::UnitRange, krange::UnitRange) = ScalarData(
#     name = data.name,
#     time = data.time,
#     grid = Grid{typeof(data.grid.scalex)}(
#         nx = irange[end]-irange[1],
#         ny = jrange[end]-jrange[1],
#         nz = krange[end]-krange[1],
#         x = data.grid.x[irange],
#         y = data.grid.y[jrange],
#         z = data.grid.z[krange],
#         lx = data.grid.x[irange[2]] - data.grid.x[irange[1]],
#         ly = data.grid.y[jrange[2]] - data.grid.x[jrange[1]],
#         lz = data.grid.z[krange[2]] - data.grid.x[krange[1]]
#     ),
#     field = view(data.field, irange, jrange, krange)
# )


end