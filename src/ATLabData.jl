"""
Module to load, process and visualise data from ATLab in Julia.
"""
module ATLabData

include("DataStructures.jl")
using .DataStructures
    export Grid, ScalarData, VectorData, AveragesData

include("IO.jl")
using .IO
    export load, Grid_from_file, search_inifile, VAR

include("Basics.jl")
using .Basics
    export size, display, +, -, *, ^, abs, log, convert, eltype
    export crop, norm, logarithm
    export component

include("Analysis.jl")
using .Analysis
    export gradient, curl

include("Statistics.jl")
using .Statistics
    export average, rms, mean, flucs, wave, turbulence

include("Visualization.jl")
using .Visualization
    export visualize
    export animate
    export heatmap, lines
    export save_for_LaTeX
    export save_visualize
    export add_BgPlot!

include("Physics.jl")
using .Physics
    export vorticity, enstrophy, Ri, tke


function __init__()
    # 
end


end