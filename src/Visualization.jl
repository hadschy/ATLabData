module Visualization

import CairoMakie: heatmap, lines

using CairoMakie
using LaTeXStrings
using ..DataStructures: 
    ScalarData, AveragesData
using ..IO: load

include("Visualization/Data2DVisualization.jl")
include("Visualization/AverageDataVisualization.jl")

export visualize
export heatmap, lines
export animate
export add_BgPlot!

end