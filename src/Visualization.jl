module Visualization

using CairoMakie
# using GLMakie

using LaTeXStrings
using ..DataStructures: 
    Data, AveragesData
using ..IO: load

include("Visualization/Data2DVisualization.jl")
include("Visualization/AverageDataVisualization.jl")

export visualize
export animate
export display_visualize
export save_visualize
export save_for_LaTeX
export add_BgPlot!

end