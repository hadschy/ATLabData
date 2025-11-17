module Tools

using Interpolations
using ..DataStructures
using ..IO: init

export shiftgrid!, transform_grid
export search_inifile

"""
Transform the grid of _data_ in _grid_. _shift_ corresponds to the axes given 
in _shiftaxis_.
"""
function transform_grid(
        data::ScalarData,
        grid::Grid;
        shift::Vector{<:AbstractFloat}=[0.0],
        shiftaxis::Vector{Symbol}=[:z]
    )
    return _transform_grid(data, grid, shift, shiftaxis)
end


function transform_grid(
        data::ScalarData,
        gridfile::String;
        shift::Vector{<:AbstractFloat}=[0.0],
        shiftaxis::Vector{Symbol}=[:z]
    )
    return _transform_grid(data, loadgrid(gridfile), shift, shiftaxis)
end


@inline function _transform_grid(
        data::ScalarData, grid::Grid, shift::Vector{<:AbstractFloat}, shiftaxis::Vector{Symbol}
    )::ScalarData
    # Container for transformed data
    newdata = init(grid)
    newdata.name = data.name
    newdata.time = data.time
    # Shift original grid
    for (i, axis) ∈ enumerate(shiftaxis)
        shiftgrid!(data, shift[i], axis=axis)
    end
    # Use interpolation of original data with higher resolution ot fill the container with the lower reolution grid
    if data.grid.ny==1
        itp = interpolate(
            (data.grid.x, data.grid.z),                 # Nodes of the grid
            data.field[:,1,:],                          # Field to be interpolated
            Gridded(Linear())                           # Interpolation type
        )
    else
        itp = interpolate(
            (data.grid.x, data.rgid.y, data.grid.z),    # Nodes of the grid
            data.field[:,:,:],                          # Field to be interpolated
            Gridded(Linear())                           # Interpolation type
        )
    end
    for k ∈ eachindex(newdata.grid.z)
        for j ∈ eachindex(newdata.grid.y)
            for i ∈ eachindex(newdata.grid.x)
                newdata.field[i,j,k] = itp(newdata.grid.x[i], newdata.grid.z[k])
            end
        end
    end
    return newdata
end


function shiftgrid!(data::ScalarData, shift::AbstractFloat; axis::Symbol=:z)
    if axis==:z
        data.grid.z .= data.grid.z .+ shift
    elseif axis==:y
        data.grid.y .= data.grid.y .+ shift
    elseif axis==:x
        data.grid.x .= data.grid.x .+ shift
    end
end


function search_inifile(file::String, block::String, key::String)::String
    f = open(file, "r")
    res = ""
    corr_block = false
    while ! eof(f)
        s = readline(f)
        # Correct block?
        if startswith(s, "[")
            if occursin(block, s)
                corr_block = true
            else
                corr_block = false
            end
        end
        # Check for key if in correct block
        if corr_block && ! startswith(s, "[")
            val = split(s, "=")
            if length(val)==2
                if key==val[1]
                    res = val[2]
                end
            end
        end
    end
    close(f)
    return res
end


end