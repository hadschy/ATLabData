module IO

using NetCDF
using ..DataStructures

export load, load! 
export file_for_time, Grid_from_file, search_inifile, VAR


"""
    load(file) -> ScalarData
Load the data contained in the path _file_ into the type _Data_.

    load(dir, field, time) -> ScalarData
Search in _dir_ for the file containing _field_ at _time_ and load into the type 
_Data_.

    load(xfile, yfile, zfile) -> VectorData
Load the data contained in the paths _xfile_, _yfile_ and _zfile_ into the 
type _VectorData_.

    load(file, var) -> AveragesData
Load the data contained in the path _file_ into the type _AveragesData_.
_file_ has to be NetCDF file containing the averages from _average.x_.
"""
load(file::String)::ScalarData = ScalarData_from_file(file)
load(dir::String, field::String, time::Real; component::String=".0")::ScalarData = load(file_for_time(dir, time, field, component))

load(xfile::String, yfile::String, zfile::String)::VectorData = VectorData_from_files(xfile, yfile, zfile)

load(file::String, var::String)::AveragesData = AveragesData_from_NetCDF(file, var)
load(dir::String, var::String, time::Real, avg::Bool) = load(avgfile_for_time(dir, time), var)

"""
    load!(data, file)
Loads the content in _file_ into preallocated _data_ of type ScalarData.
Assumes that _data_containes an already loaded _grid_.
"""
function load!(data::ScalarData, file::String)
    ScalarData_from_file!(data, file)
end


# ------------------------------------------------------------------------------
#                                   Grid
# ------------------------------------------------------------------------------
""" 
    Load the file _grid_ into composite type _Grid_
"""
function Grid_from_gridfile(gridfile::String)::Grid
    "f90_record_markers from the amount of write commands in storing routine "
    io = open(gridfile, "r")
    marker1 = read(io, Int32)
    nx = read(io, Int32)
    ny = read(io, Int32)
    nz = read(io, Int32)
    marker2 = read(io, Int32)
    marker3 = read(io, Int32)
    scalex = read(io, Float64)
    scaley = read(io, Float64)
    scalez = read(io, Float64)
    marker4 = read(io, Int32)
    marker5 = read(io, Int32)
    x = Vector{Float64}(undef, nx)
    for i ∈ 1:nx
        x[i] = read(io, Float64)
    end
    marker7 = read(io, Int32)
    marker8 = read(io, Int32)
    y = Vector{Float64}(undef, ny)
    for i ∈ 1:ny
        y[i] = read(io, Float64)
    end
    marker9 = read(io, Int32)
    marker10 = read(io, Int32)
    z = Vector{Float64}(undef, nz)
    for i ∈ 1:nz
        z[i] = read(io, Float64)
    end
    return Grid(nx, ny, nz, scalex, scaley, scalez, x, y, z)
end


function Grid_from_file(dir::String)::Grid
    gridfile = joinpath(dir, "grid")
    return Grid_from_gridfile(gridfile)
end


# ------------------------------------------------------------------------------
#                               ScalarData
# ------------------------------------------------------------------------------
"""
    Loading data that is stored in a single file.
"""
function ScalarData_from_file(fieldfile::String)::ScalarData
    verbose("ScalarData", fieldfile)
    filename = split(fieldfile, "/")[end]
    if startswith(filename, "flow.") || startswith(filename, "scal.")
        return ScalarData_from_raw(fieldfile)
    else
        return ScalarData_from_visuals(fieldfile)
    end
end


function ScalarData_from_raw(fieldfile::String)::ScalarData{Float64, Int32}
    grid = Grid_from_file(dirname(fieldfile))
    buffer, t = Array_from_rawfile(grid, fieldfile)
    return ScalarData(
        name = splitpath(fieldfile)[end],
        grid = grid,
        time = t,
        field = buffer
    )
end


function ScalarData_from_visuals(fieldfile::String)::ScalarData{Float32, Int32}
    grid = convert(Float32, Grid_from_file(dirname(fieldfile)))
    return ScalarData(
        name = splitpath(fieldfile)[end],
        grid = grid,
        time = time_from_file(fieldfile),
        field = Array_from_file(grid, fieldfile)
    )
end


function ScalarData_from_file!(data::ScalarData, fieldfile::String)
    verbose("ScalarData", fieldfile)
    filename = split(fieldfile, "/")[end]
    if startswith(filename, "flow.") || startswith(filename, "scal.")
        (data.field, data.time) = Array_from_rawfile(data.grid, fieldfile)
    else
        data.time = time_from_file(fieldfile)
        data.field .= Array_from_file(data.grid, fieldfile)
    end
end


# ------------------------------------------------------------------------------
#                               VectorData
# ------------------------------------------------------------------------------
"""
    Loading data from three files. No check for physical consistency provided.
"""
function VectorData_from_files(
        xfieldfile::String,
        yfieldfile::String,
        zfieldfile::String
    )::VectorData
    verbose("VectorData", xfieldfile, yfieldfile, zfieldfile)
    filename = split(xfieldfile, "/")[end]
    if startswith(filename, "flow.") || startswith(filename, "scal.")
        return VectorData_from_raw(xfieldfile, yfieldfile, zfieldfile)
    else
        return VectorData_from_visuals(xfieldfile, yfieldfile, zfieldfile)
    end
end


function VectorData_from_raw(
        xfieldfile::String,
        yfieldfile::String,
        zfieldfile::String
    )::VectorData{Float64, Int32}
    grid = Grid_from_file(dirname(xfieldfile))
    buffer, t = Array_from_rawfile(grid, xfieldfile)
    return VectorData(
        name = string(splitpath(xfieldfile)[end][1:end-2]),
        grid = grid,
        time = t,
        xfield = buffer,
        yfield = Array_from_rawfile(grid, yfieldfile)[1],
        zfield = Array_from_rawfile(grid, zfieldfile)[1]
    )
end


function VectorData_from_visuals(
        xfieldfile::String,
        yfieldfile::String,
        zfieldfile::String
    )
    grid = convert(Float32, Grid_from_file(dirname(xfieldfile)))
    return VectorData(
        name = string(splitpath(xfieldfile)[end][1:end-2]),
        grid = grid,
        time = time_from_file(xfieldfile),
        xfield = Array_from_file(grid, xfieldfile),
        yfield = Array_from_file(grid, yfieldfile),
        zfield = Array_from_file(grid, zfieldfile)
    )
end


# ------------------------------------------------------------------------------
#                         Array from binary
# ------------------------------------------------------------------------------
"""
    Load array from a binary file according to the information in _grid_.
"""
function Array_from_file(
        grid::Grid{T,I}, fieldfile::String
    )::Array{T,3} where {T<:AbstractFloat, I<:Signed}
    buffer = Vector{T}(undef, grid.nx*grid.ny*grid.nz)
    read!(fieldfile, buffer)
    return reshape(buffer, (grid.nx, grid.ny, grid.nz))
end
 

function Array_from_rawfile(
        grid::Grid{T,I}, fieldfile::String
    )::Tuple{Array{T, 3}, Float64} where {T<:AbstractFloat, I<:Signed}
    io = open(fieldfile, "r")
    """ 
        Header contains: 
        headersize/offset, nx, ny, nz, nt, time, visc, froude, schmidt
        unformatted with no record markers
    """
    headersize = read(io, Int32) # First entry is the header length in bytes
    seek(io, 5*sizeof(headersize))
    time = read(io, Float64)
    seek(io, headersize) # Jump to first last entry belonging to the header
    buffer = Vector{Float64}(undef, grid.nx*grid.ny*grid.nz)
    read!(io, buffer)
    return reshape(buffer, (grid.nx, grid.ny, grid.nz)), time
end


# ------------------------------------------------------------------------------
#                            AveragesData
# ------------------------------------------------------------------------------
function VAR(file::String, time::Real, var::String)::Vector{Float32}
    times = ncread(file, "t")
    Δt = time
    index = 1
    for i ∈ eachindex(times)
        if abs(times[i] - time) < Δt
            index = i
            Δt = abs(times[i] - time)
        end
    end
    return ncread(file, var)[:,index]
end


function AveragesData_from_allNetCDF(;
        file::String, 
        var::String, 
        time::Float64
    )::AveragesData
    verbose("AveragesData", file)
    t = ncread(file, "t")
    range = ncread(file, "y")
    data = ncread(file, var)
    # Find index for time in t:
    minΔt = 1000000000.
    index = 0
    for i in eachindex(t)
        if abs(t[i]-time) < minΔt
            index = i
            minΔt = abs(t[i]-time)
        end
    end
    println("Found averages data for: t = ", t[index])
    println("Was lokking for: time = ", time)
    return AveragesData(name=var, time=time, field=data[:,index], range=range)
end


function AveragesData_from_NetCDF(
        file::String, 
        var::String
    )::AveragesData
    verbose("AveragesData", file)
    t = ncread(file, "t")
    range = ncread(file, "z")
    data = ncread(file, var)
    return AveragesData(name=var, time=t[1], field=data[:,1], range=range)
end


function avgfile_for_time(dir::String, time::Real)::String
    println("Looking for t = $(time) ...")
    filenames = filter(x -> startswith(x, "avg"), readdir(dir, join=false))
    filenames = filter(x -> !occursin("_all", x), filenames)
    Δt = time
    Δt_max = 0.1*time
    return_file = ""
    for filename ∈ filenames
        file = joinpath(dir, filename)
        t = ncread(file, "t")[1]
        δt = abs(t - time)
        if δt < Δt
            return_file = file
            Δt = δt
        end
    end
    if Δt > Δt_max
        @warn "Maximal relative error for _time_ exceeded!"
    end
    if return_file==""
        error("No avg file found for in $(dir).")
    end
    printstyled("   $(return_file)     $(ncread(return_file, "t")[1]) \n", color=:cyan)
    return return_file
end


function AveragesData_to_txtfile(data::AveragesData, outdir::String)
    file = joinpath(
        outdir, "avgs_"*string(data.name)*"_"*string(data.time)*".txt"
    )
    verbose("AveragesData", file, save=true)
    f = open(file, "w")
    for i in eachindex(data.range)
        println(f, data.field[i], "    ", data.range[i])
    end
    close(f)
    return nothing
end


# ------------------------------------------------------------------------------
#                             Helping functions
# ------------------------------------------------------------------------------
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


function timestep_from_filename(filename::String)::Int
    namestring = split(filename, "/")[end]
    namestring = split(namestring, ".")[1]
    stepstring = namestring[end-5:end]
    return parse(Int, stepstring)
end


function time_from_timestep(timestep::Int, dir::String)::AbstractFloat
    avgfile = "avg"*string(timestep)*".nc"
    avgfile = joinpath(dir, avgfile)
    return ncread(avgfile, "t")[1]
end


function time_from_file(file::String)::Float32
    tstep = timestep_from_filename(file)
    t = time_from_timestep(tstep, dirname(file))
    return t
end


function file_for_time(
        dir::String, 
        time::Real,
        field::String,
        component::String=".0"
    )::String
    Δt = time
    f = ""
    maxΔt_relative = 0.1
    maxΔt = maxΔt_relative*time
    println("Searching for t = ", time)
    for filename in readdir(dir)
        correct_field = false
        if component ∈ (".1", ".2", ".3", ".4", ".5", ".6")
            if startswith(filename, field) && endswith(filename, component)
                correct_field = true
            end
        elseif component == ".0"
            if startswith(filename, field)
                correct_field = true
            end
        else
            error("Loading a field with this _component_ is not implemented.")
        end
        if correct_field
            file = joinpath(dir, filename)
            t = time_from_file(file)
            printstyled("   ", filename, "    ", t, "\n", color=:cyan)
            if abs(time - t) < Δt
                f = file
                Δt = abs(time - t)
            end
        end
    end
    if Δt > maxΔt
        @warn "The time loaded exceeds the maximal accepted relative error of $(maxΔt_relative)."
    end
    return f
end


# """ 
# Search for the proper timestep that is nearest to _time_ based on _avg_all.nc_ 
# (from ncrcat). Return the correponding field string.
# """
# function file_for_time_new(
#         dir::String, 
#         time::Float64, 
#         field::String, 
#         component::String=".0"
#     )::String
#     println("Searching for t = $(time)")
#     t, avgfile = time_from_avgall(dir, time)
#     timestep = split(avgfile, "avg", keepempty=false)[1]
#     timestep = split(timestep, ".")[1]
#     for i ∈ 1:(6-length(timestep))
#         timestep = "0"*timestep
#     end
#     filename = field * timestep
#     if component ∈ (".1", ".2", ".3")
#         filename = filename * component
#     end
#     return filename
# end


# function time_from_avgall(dir::String, time::Float64)::Tuple{Float64, String}
#     Δt = time
#     t = ncread(joinpath(dir, "avg_all.nc"), "t")
#     itime = 0
#     for i in eachindex(t)
#         if abs(time-t[i]) < Δt
#             Δt = abs(time-t[i])
#             itime = i
#         end
#     end
#     avgfile = 
#     return t[itime], avgfile
# end


function verbose(dtype::String, files...; save::Bool=false)
    if !save
        text = "Loading " * dtype * " ..."
    else
        text = "Saving " * dtype * " ..."
    end
    println(text)
    for file in files
        printstyled("   "*file*"\n", color=:cyan)
    end
end


function done()
    printstyled("   Done \n", color=:green)
end


end