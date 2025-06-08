"""
Module to load, process and visualise data from ATLab in Julia.
"""
module ATLabData

include("DataStructures.jl")
using .DataStructures
    export Grid, ScalarData, VectorData, AveragesData

include("IO.jl")
using .IO
    export load, Grid_from_file, search_inifile

include("DataOperations.jl")
using .DataOperations
    export size, display, +, -, *, ^, abs, log, convert, eltype
    export crop, gradient, norm, curl, logarithm
    export component
    export rms, average

include("Visualization.jl")
using .Visualization
    export visualize
    export animate
    export display_visualize
    export save_for_LaTeX
    export save_visualize
    export add_BgPlot!

include("Physics.jl")
using .Physics
    export vorticity, enstrophy, Ri


function __init__()
    # 
end


end