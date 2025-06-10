function lines(
        data::AveragesData;
        xlabel::String=data.name,
        ylabel::String="z",    
    )::Tuple{Figure, Axis, Lines}
    fig = Figure(size=(1000, 1000))
    ax = Axis(
        fig[1,1],
        titlesize = 25,
        xlabel = xlabel, ylabel = ylabel, 
        ylabelsize = 40, xlabelsize = 40, 
        xgridstyle = :dash, 
        ygridstyle = :dash, 
        xtickalign = 1, xticksize = 10, 
        ytickalign = 1, yticksize = 10,
        xticklabelsize = 40, yticklabelsize = 40, 
        xlabelpadding = 10, ylabelpadding = 10
    )
    ln = lines!(ax, view(data.field, :), view(data.grid.z, :), label="t = "*string(data.time))
    axislegend(ax)
    return fig, ax, ln
end


function visualize(
        data::AveragesData; 
        xlabel = "$(data.name)",
        ylabel = "z",
    )
    fig, ax, ln = lines(data; xlabel=xlabel, ylabel=ylabel)
    display(fig)
end