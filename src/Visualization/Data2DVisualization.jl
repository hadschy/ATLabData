# TODO: proper resolution scaling


function save_for_LaTeX(
        data::Data, 
        outdir::String; 
        slice::Int=1,
        sizey::Int=1000
    )
    
    # sizez = round(Int32, data.grid.nx/data.grid.nz*sizezx

    set_theme!(theme_latexfonts())
    fig = Figure(
    size=(sizex, sizez),
    figure_padding=0., 
    # backgroundcolor=:transparent
    )
    ax = Axis(
        fig[1,1], 
        xlabelsize=30, ylabelsize=30,
        xticklabelsize=25, yticklabelsize=25,
        xlabel=L"x", ylabel=L"z",
        # title=L"u_y",
        titlesize=40
    )
    # Note: PGFplots only knows few colormaps, therefore using viridis here
    heatmap!(
        ax,
        data.grid.x,
        data.grid.z,
        data.field[:,slice,:], 
        colormap=:viridis
    )
    hidedecorations!(ax)
    outfile = joinpath(outdir, data.name*"-pgfplots.png")
    save(outfile, fig)

    # Save time information
    t = data.time
    println(t)
    filename = data.name*".time"
    outfile = joinpath(outdir, filename)
    println("Storing physical time:    ", outfile)
    f = open(outfile, "w")
    write(f, string(t)*"\n")
    close(f)

    # Save data range (e.g. for colorbar)
    filename = data.name * ".range"
    outfile = joinpath(outdir, filename)
    println("Storing data range:    ", outfile)
    f = open(outfile, "w")
    write(f, string(maximum(data.field)) * "\n")
    write(f, string(minimum(data.field)) * "\n")
    close(f)
end


function save_for_LaTeX(fig::Figure, outfile::String)
    fig.figure_padding = 0
    hidedecorations!(ax)
    save(fig, outfile)
end


"""
    visualize(
        data; 
        slice=1, 
        sizex=1000, 
        colormap=:diff, 
        colorrange_max=maximum(data.field),
        colorrange_min=minimum(data.field),
        label="default",
        fontsizescaling=1,
        critical_level=-1.0,
        xlabel="x", ylabel="z",
        xtickstep=10, ytickstep=20,
        scaling=1.0, cscale=identity,
        title=data.name*"  ;  t = "*string(data.time)
    ) -> Figure, Axis
Main visualization routine, that generates and returns the figure and axis by 
utilizing metadata of the _Data_ type. The returned axis contains a heatmap 
of _data.field[:,:,slice]_.
"""
function visualize(
        data::Data;
        slice::Int=1, 
        sizex::Int=2000,
        colormap=:RdGy,
        colorrange_max::Float32=maximum(data.field),
        colorrange_min::Float32=minimum(data.field),
        label=" ",
        fontsize::Int=round(Int, sizex*40/2000),
        fontsizescaling=1,
        critical_level::Float64=-1.0,
        xlabel=L"x",
        ylabel=L"z",
        xtickstep=10, #round(data.grid.scalex/8),
        ytickstep=10, #round(data.grid.scaley/4),
        scaling::Float64=1.0,
        cscale=identity,
        title::String=data.name*"  ;  t = "*string(data.time),
        colorbarticks::Vector{Float32}=[
            round(colorrange_min, digits=2), 
            round((colorrange_max+colorrange_min)/2, digits=2), 
            round(colorrange_max, digits=2)
        ],
        # interactive=true
    )::Tuple{Figure, Axis, Heatmap}
    println("Visualizing ...")
    printstyled("   $(data.name) \n", color=:cyan)
    sizez = round(Int32, data.grid.scalez/data.grid.scalex*sizex)    
    # if interactive
    #     printstyled("   Backend: GLMakie \n", color=:gray, italic=true)
    #     GLMakie.activate!()
    #     fig = Figure(size=(sizex, sizez))
    #     ax = Axis(
    #         fig[1,1], 
    #         xlabel=xlabel, ylabel=ylabel,
    #         title=title,
    #         backgroundcolor=:transparent
    #     )
    # else
        printstyled("   Backend: CairoMakie \n", color=:gray, italic=true)
        CairoMakie.activate!()
        fontsize = round(Int, fontsizescaling*fontsize)
        fig = Figure(
            size=(sizex, sizez),
            px_per_iunit=10,
            # backgroundcolor=:transparent,
            figure_padding=round(Int, 10/2000*sizex)
        )
        ax = Axis(
            fig[1,1], 
            xlabel=xlabel, ylabel=ylabel,
            xlabelsize=fontsize, ylabelsize=fontsize,
            xticklabelsize=fontsize, yticklabelsize=fontsize,
            xticks=data.grid.x[1]:xtickstep:data.grid.x[end],
            yticks=round(data.grid.z[1]):ytickstep:round(data.grid.z[end]),
            title=title,
            titlesize=fontsize/2,
            backgroundcolor=:transparent
        )
        rowsize!(fig.layout, 1, Aspect(1, data.grid.nz/data.grid.nx))
        resize_to_layout!(fig)
    # end
    hm = heatmap!(
        ax,
        view(data.grid.x ./ scaling, :), 
        view(data.grid.z ./ scaling, :),
        view(data.field, :, slice, :),
        colormap=colormap,
        colorrange=(colorrange_min, colorrange_max),
        colorscale=cscale
    )
    # if interactive
    #     Colorbar(fig[1,2], hm, label="")
    # else
        Colorbar(
            fig[1,2], 
            hm, 
            ticklabelsize=Int(fontsize), 
            label=label,
            labelsize=Int(fontsize),
            size = round(Int, sizex*12/2000),
            ticks=colorbarticks
        )
    # end
    if critical_level > 0.0
        lines!(
            ax,
            view(data.grid.x ./ scaling, :),
            ones(data.grid.nx)*critical_level
        )
    end
    return fig, ax, hm
end


"""
    add_BgPlot!(
        fig, data;
        xlabel=L"u_0",
        fontsizescaling=1,
        xmin=minimum(data.field),
        mxmax=maximum(data.field)
    ) -> Figure
Adds a line plot of data to fig[1,2]. Intented to be used with ```visualize````
"""
function add_BgPlot!(
        fig::Figure, 
        data::AveragesData;
        xlabel::LaTeXString=L"u_0",
        fontsizescaling=1,
        xmin=minimum(data.field),
        xmax=maximum(data.field),
        xticks=[round(minimum(data.field), digits=1), round(maximum(data.field), digits=1)],
        yticks=[round(data.range[1]), round(data.range[end])]
    )::Figure
    width = round(Int, fig.scene.viewport[].widths[1] * 100/2000)
    fontsize = round(Int, width * 40/100)
    fontsize = round(Int, fontsizescaling*fontsize)
    linewidth = round(Int, width * 5/100)
    ax = Axis(
        fig[1,3];
        xlabel=xlabel, ylabel= L"z",
        xlabelsize=fontsize, ylabelsize=fontsize,
        xticklabelsize=fontsize, yticklabelsize=fontsize,
        xticks=xticks,
        yticks=yticks,
        limits=(xmin, xmax, data.range[1], data.range[end]),
        backgroundcolor=:transparent,
        width=width
    )
    lines!(ax, data.field, data.range, linewidth=linewidth)
    # hideydecorations!(ax)
    resize_to_layout!(fig)
    return fig
end


function display_visualize(data::Data; slice::Int=1)
    set_theme!(merge(theme_dark(), theme_latexfonts()))
    fig, ax, hm = visualize(data, slice=slice)
    display(fig)
end


function animate(
        dir::String, 
        field::String;
        fps::Int=2,
        loader::Function=load,
        visualizer::Function=visualize,
        live::Bool=false,
        outfile::String=joinpath(dir, "video/$(field).mp4"),
        slice::Int=1
    )
    println("Animating:")
    printstyled("   $(dir)", color=:cyan)
    filenames = filter(x -> startswith(x, field), readdir(dir, join=false))
    file = joinpath(dir, filenames[1])
    data = loader(file)
    fig, ax, hm = visualizer(data)
    # if live
    #     # Live in GLWindow using GLMakie
    #     display(fig)
    #     for i ∈ eachindex(filenames)
    #         file = joinpath(dir, filenames[i])
    #         data = loader(file)
    #         hm[3] = view(data.field, :, slice, :)
    #         ax.title = "$(field) ; t = $(data.time)"
    #         sleep(1/fps)
    #     end
    # else
        # Save as .mp4 file with CairoMakie
        if ! ispath(dirname(outfile))
            mkpath(dirname(outfile))
        end
        record(fig, outfile, 1:length(filenames), framerate=fps) do i
            file = joinpath(dir, filenames[i])
            data = loader(file)
            hm[3] = view(data.field, :, slice, :)
            ax.title = "$(field) ; t = $(data.time)"
        end
    # end
end