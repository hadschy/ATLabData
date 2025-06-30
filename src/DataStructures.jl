module DataStructures

export AbstractData
export Grid, ScalarData, VectorData, AveragesData


""" Parent abstract type for all composite types containing the data """
abstract type AbstractData{T<:AbstractFloat, I<:Signed} end


struct Grid{T,I} <: AbstractData{T,I}
    nx::I
    ny::I
    nz::I
    scalex::T
    scaley::T
    scalez::T
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
end
Grid(;
    nx::Signed, ny::Signed, nz::Signed, 
    lx::AbstractFloat, ly::AbstractFloat, lz::AbstractFloat,
    x::Vector{<:AbstractFloat}, 
    y::Vector{<:AbstractFloat}, 
    z::Vector{<:AbstractFloat}
) = Grid(nx, ny, nz, lx, ly, lz, x, y, z)


struct ScalarData{T,I} <: AbstractData{T,I}
    name::String
    grid::Grid{T,I}
    time::T
    field::Array{T,3}
end
ScalarData(;
    name::String, grid::Grid, time::AbstractFloat, field::Array{<:AbstractFloat, 3}
) = ScalarData(name, grid, time, field)


struct VectorData{T,I} <: AbstractData{T,I}
    name::String
    grid::Grid{T,I}
    time::T
    xfield::Array{T,3}
    yfield::Array{T,3}
    zfield::Array{T,3}
end
VectorData(;
    name::String, grid::Grid, time::AbstractFloat, 
    xfield::Array{<:AbstractFloat,3}, 
    yfield::Array{<:AbstractFloat,3}, 
    zfield::Array{<:AbstractFloat,3},
) = VectorData(name, grid, time, xfield, yfield, zfield)


struct AveragesData{T,I} <: AbstractData{T,I}
    name::String
    time::T
    # range::Vector{T}
    grid::Grid{T,I}
    field::Vector{T}
end
AveragesData(;
    name::String,
    time::AbstractFloat,
    range::Vector{<:AbstractFloat},
    field::Vector{<:AbstractFloat}
) = AveragesData(
        name, time, 
        Grid{eltype(time), Int}(
            1, 1, length(range), 
            0.0, 0.0, range[end],
            [0.0], [0.0], range
        ),
        field
    )


function Base.getindex(data::ScalarData, i::Int, j::Int, k::Int)
    return data.field[i,j,k]
end


#-------------------------------------------------------------------------------
#                            Under  Construction
#-------------------------------------------------------------------------------
# TODO SubData in order ot realise perfomance with view(data)



# """ Needed for defining a view method for Data """
# struct SubData{T,I} <: AbstractData{T,I}
#     name::String,
#     time::T
#     grid::Grid{T,I},
#     field::SubArray
# end


end