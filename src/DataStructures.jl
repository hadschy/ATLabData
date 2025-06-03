module DataStructures

export Grid, Data, ScalarData, VectorData, AveragesData



# struct Grid{T<:AbstractFloat}
#     nx::Int32
#     ny::Int32
#     nz::Int32
#     scalex::T
#     scaley::T
#     scalez::T
#     dx::T
#     dy::T
#     dz::T
#     x::Vector{T}
#     y::Vector{T}
#     z::Vector{T}
#     function Grid{T}(;
#             nx::Int32, ny::Int32, nz::Int32, 
#             lx::T, ly::T, lz::T, 
#             x::Vector{T}, y::Vector{T}, z::Vector{T}
#         ) where {T<:AbstractFloat}
#         return new(nx, ny, nz, lx, ly, lz, lx/nx, ly/ny, lz/nz, x, y, z)
#     end
# end


# struct Data
#     name::String
#     grid::Grid{Float64}
#     field::Array{Float32, 3}
#     time::Float32
#     function Data(;
#             grid::Grid{Float64}, 
#             field::Array{Float32, 3},
#             name::String="DefaultName",
#             time::Float32
#         )
#         return new(name, grid, field, time)
#     end
# end


# struct VectorData
#     name::String
#     grid::Grid
#     time::Float32
#     xfield::Array{Float32, 3}
#     yfield::Array{Float32, 3}
#     zfield::Array{Float32, 3}
#     function VectorData(;
#             name::String="DefaultName",
#             grid::Grid,
#             time::Float32,
#             xfield::Array{Float32, 3},
#             yfield::Array{Float32, 3},
#             zfield::Array{Float32, 3}
#         )
#         return new(name, grid, time, xfield, yfield, zfield)
#     end
# end


# struct AveragesData
#     name::String
#     time::Float32
#     range::Vector{Float32}
#     field::Vector{Float32}
#     function AveragesData(;
#             name::String="DefaultName",
#             time::Float32,
#             field::Vector{Float32},
#             range::Vector{Float32}
#         )
#         return new(name, time, range, field)
#     end
# end

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
    function Grid{T, I}(;
            nx::I, ny::I, nz::I, lx::T, ly::T, lz::T, 
            x::Vector{T}, y::Vector{T}, z::Vector{T}
        ) where {T<:AbstractFloat, I<:Signed}
        return new(nx, ny, nz, lx, ly, lz, x, y, z)
    end
end


struct Data
    name::String
    grid::Grid{Float64, Int32}
    field::Array{Float32, 3}
    time::Float32
    function Data(;
            grid::Grid{Float64, Int32}, 
            field::Array{Float32, 3},
            name::String="DefaultName",
            time::Float32
        )
        return new(name, grid, field, time)
    end
end


struct ScalarData{G<:AbstractFloat,T,I} <: AbstractData{T,I}
    name::String
    grid::Grid{G,I}
    time::T
    field::Array{T,3}
    function ScalarData{G,T,I}(;
            name::String, time::T, grid::Grid{G,I}, field::Array{T,3}
        ) where {G<:AbstractFloat, T<:AbstractFloat, I<:Signed}
        return new(name, grid, time, field)
    end
end


struct VectorData{T,I} <: AbstractData{T,I}
    name::String
    grid::Grid{T,I}
    time::T
    xfield::Array{T,3}
    yfield::Array{T,3}
    zfield::Array{T,3}
    function VectorData{T,I}(;
            name::String="DefaultName",
            grid::Grid{T,I},
            time::T,
            xfield::Array{T,3},
            yfield::Array{T,3},
            zfield::Array{T,3}
        ) where {T<:AbstractFloat, I<:Signed}
        return new(name, grid, time, xfield, yfield, zfield)
    end
end


struct AveragesData{T,I} <: AbstractData{T,I}
    name::String
    time::T
    range::Vector{T}
    field::Vector{T}
    function AveragesData{T,I}(;
            name::String="DefaultName",
            time::T,
            field::Vector{T},
            range::Vector{T}
        ) where {T<:AbstractFloat, I<:Signed}
        return new(name, time, range, field)
    end
end



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