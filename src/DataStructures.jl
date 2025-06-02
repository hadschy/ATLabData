module DataStructures

export Grid, Data, VectorData, AveragesData



struct Grid{T<:AbstractFloat}
    nx::Int32
    ny::Int32
    nz::Int32
    scalex::T
    scaley::T
    scalez::T
    dx::T
    dy::T
    dz::T
    x::Vector{T}
    y::Vector{T}
    z::Vector{T}
    function Grid{T}(;
            nx::Int32, ny::Int32, nz::Int32, 
            lx::T, ly::T, lz::T, 
            x::Vector{T}=collect(range(0.0, lx, nx)), 
            y::Vector{T}=collect(range(0.0, ly, ny)), 
            z::Vector{T}=collect(range(0.0, lz-lz/nz, nz))
        ) where {T<:AbstractFloat}
        return new(nx, ny, nz, lx, ly, lz, lx/nx, ly/ny, lz/nz, x, y, z)
    end
end


struct Data
    name::String
    grid::Grid{Float64}
    field::Array{Float32, 3}
    time::Float32
    function Data(;
            grid::Grid{Float64}, 
            field::Array{Float32, 3},
            name::String="DefaultName",
            time::Float32
        )
        return new(name, grid, field, time)
    end
end


struct VectorData
    name::String
    grid::Grid
    time::Float32
    xfield::Array{Float32, 3}
    yfield::Array{Float32, 3}
    zfield::Array{Float32, 3}
    function VectorData(;
            name::String="DefaultName",
            grid::Grid,
            time::Float32,
            xfield::Array{Float32, 3},
            yfield::Array{Float32, 3},
            zfield::Array{Float32, 3}
        )
        return new(name, grid, time, xfield, yfield, zfield)
    end
end


struct AveragesData
    name::String
    time::Float32
    range::Vector{Float32}
    field::Vector{Float32}
    function AveragesData(;
            name::String="DefaultName",
            time::Float32,
            field::Vector{Float32},
            range::Vector{Float32}
        )
        return new(name, time, range, field)
    end
end


function Base.getindex(data::Data, i::Int, j::Int, k::Int)
    return data.field[i,j,k]
end


#-------------------------------------------------------------------------------
#                            Under  Construction
#-------------------------------------------------------------------------------
# TODO AbstractData and SubData in order ot realise perfomance with view(data)

# """ Parent abstract type for all composite types containing the data """
# abstract type AbstractData{T<:AbstractFloat, I<:Int} end


# # To be the replacment for Grid
# struct NewGrid{FloatT, InT} <: AbstractData{FloatT, IntT}
#     nx::IntT
#     ny::IntT
#     nz::IntT
#     scalex::FloatT
#     scaley::FloatT
#     scalez::FloatT
#     x::Vector{FloatT}
#     y::Vector{FloatT}
#     z::Vector{FloatT}
#     function NewGrid{FloatT, IntT}(;
#             nx::IntT, ny::IntT, nz::IntT, 
#             lx::FloatT, ly::FloatT, lz::FloatT, 
#             x::Vector{FloatT}=collect(FloatT, range(0.0, lx, nx)), 
#             y::Vector{FloatT}=collect(FloatT, range(0.0, ly, ny)), 
#             z::Vector{FloatT}=collect(FloatT, range(0.0, lz-lz/nz, nz))
#         ) where {FloatT<:AbstractFloat, IntT<:Int}
#         return new(nx, ny, nz, lx, ly, lz, x, y, z)
#     end
# end


# # To be the replacment for Data
# struct ScalarData{T,I} <: AbstractData{T,I}
#     name::String
#     grid::NewGrid{T,I}
#     time::T
#     field::Array{T,3}
#     function NewData{T,I}(;
#             name::String, time::T, grid::NewGrid{T}, field::Array{T,3}
#         ) where {T<:AbstractFloat, I<:Int}
#         return new(name, grid, time, field)
#     end
# end


# # To be the replacement for VectroData
# struct NewVectorData{T,I} <: AbstractData{T,I}
#     name::String
#     grid::NewGrid{T,I}
#     time::T
#     xfield::Array{T,3}
#     yfield::Array{T,3}
#     zfield::Array{T,3}
#     function NewVectorData{T,I}(;
#             name::String="DefaultName",
#             grid::NewGrid{T,I},
#             time::T,
#             xfield::Array{T,3},
#             yfield::Array{T,3},
#             zfield::Array{T,3}
#         ) where {T<:AbstractFloat, I<:Int}
#         return new(name, grid, time, xfield, yfield, zfield)
#     end
# end


# # To be replacement for AveragesData
# struct NewAveragesData{T} <: AbstractData{T,I}
#     name::String
#     time::T
#     range::Vector{T}
#     field::Vector{T}
#     function NewAveragesData{T}(;
#             name::String="DefaultName",
#             time::T,
#             field::Vector{T},
#             range::Vector{T}
#         ) where {T<:AbstractFloat}
#         return new(name, time, range, field)
#     end
# end


# """ Needed for defining a view method for Data """
# struct SubData{T,I} <: AbstractData{T,I}
#     name::String,
#     time::T
#     grid::NewGrid{T,I},
#     field::SubArray
# end


end