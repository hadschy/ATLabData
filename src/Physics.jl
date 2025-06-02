module Physics

using ..DataStructures
using ..DataOperations
using ..IO

export vorticity, enstrophy, Ri


"""
    vorticity(data) -> VectorData
Return the curl of data, thus is a physical alternatice to curl if data 
is a velocity field.

    vorticity(dir, time) -> VectorData
Looks for the proper velocity files in dir that are nearest to _time_ and 
computes the curl.
"""
vorticity(u::VectorData)::VectorData = curl(u)
vorticity(dir::String, time::Real)::VectorData = curl(load(
    file_for_time(dir, "VelocityVector", time, ".1"),
    file_for_time(dir, "VelocityVector", time, ".2"),
    file_for_time(dir, "VelocityVector", time, ".3")
))


"""
    enstrophy(u) -> Data
Calculates the enstrophy of the givem velocity field.

    enstrophy(dir, time) -> Data
Looks in _dir_ for the velocity field at _time_ and calculates the appropriate 
    enstrophy.
"""
enstrophy(u::VectorData)::Data = Data(
    name = "enstrophy(" * u.name * ")", 
    time = u.time, 
    grid = u.grid, 
    field = norm(vorticity(u)).field.^2
)
enstrophy(dir::String, time::Real)::Data = enstrophy(load(
    file_for_time(dir, "VelocityVector", time, ".1"),
    file_for_time(dir, "VelocityVector", time, ".2"),
    file_for_time(dir, "VelocityVector", time, ".3")
))


"""
Computes and returns the local Richardson number field from the given 
buoyancy and velocity fields.

    Ri(b, u) -> Data
_u_ is given as _VectorData_.

    RI(b, ux, uy, uz) -> Data
The single components are given as _Data_.

    Ri(dir, time) -> Data
Looks in _dir_ for the buoyancy and velocity fields for _time_.
"""
Ri(b::Data, u::VectorData)::Data = Data(
    name = "Rig("*u.name*")",
    time = b.time,
    grid = b.grid,
    field = norm(gradient(u)).field.^2 ./ norm(gradient(b))
)
Ri(b::Data, ux::Data, uy::Data, uz::Data)::Data = Data(
    name = "Rig("*b.name*")",
    time = b.time,
    grid = b.grid,
    field = norm(gradient(b)).field ./ (norm(gradient(ux)).field.^2 + norm(gradient(uy)).field.^2 + norm(gradient(uz)).field.^2)
)
Ri(dir::String, time::Real)::Data = Ri(
    load(dir, "Buoyancy", time),
    load(dir, "VelocityVector", time, ".1"),
    load(dir, "VelocityVector", time, ".2"),
    load(dir, "VelocityVector", time, ".3"),
)


function Reynolds_stress(u::VectorData)::Matrix
    # TODO
end


end