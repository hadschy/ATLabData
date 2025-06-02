function visualize(
        data::AveragesData;
        xlabel::String=data.name,
        ylabel::String="z",    
    )::Tuple{Figure, Axis, Lines}
    GLMakie.activate!()
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
    ln = lines!(ax, data.field, data.range, label="t = "*string(data.time))
    axislegend(ax)
    return fig, ax, ln
end


function display_visualize(
        data::AveragesData; 
        xlabel = "$(data.name)",
        ylabel = "z",
    )
    fig, ax, ln = visualize(data; xlabel=xlabel, ylabel=ylabel)
    display(fig)
end