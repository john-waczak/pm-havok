using CairoMakie
using Dates
using CSV, DataFrames
using LinearAlgebra

include("./makie-defaults.jl")
set_theme!(mints_theme)

include("./ips7100.jl")

data_path = "./data"
bam_path = joinpath(data_path, "bam.csv")
ips_path = joinpath(data_path, "ips7100-hc.csv")

figpath = "./figures/1__exploratory-data-analysis"

if !ispath(figpath)
    mkpath(figpath)
end



# Read data
df = CSV.read(ips_path, DataFrame);
df_bam = CSV.read(bam_path, DataFrame);
dropmissing!(df)
dropmissing!(df_bam)

# https://pierasystems.com/wp-content/uploads/2024/02/IPS-Datasheet-V1.3.8.pdf
# Piera uses GRIMM EDM180 for individual device calibration for all output channels independently

# get relative humidity

# Parse dateTime column
date_fmt = dateformat"yyyy-mm-dd HH:MM:SS"
df.dateTime .= DateTime.(String.(df.dateTime), date_fmt);
df.seconds = [0.0, [ d.value ./ 1000 for d in df.dateTime[2:end] .- df.dateTime[1]]...]
df.minutes = df.seconds ./ 60.0
df.days = df.minutes ./ (60*24)


df_bam.dateTime .= DateTime.(String.(df_bam.dateTime), date_fmt);
df_bam.seconds = [0.0, [ d.value ./ 1000 for d in df_bam.dateTime[2:end] .- df_bam.dateTime[1]]...]
df_bam.minutes = df_bam.seconds ./ 60.0
df_bam.days = df_bam.minutes ./ (60*24)



ndays = (df.dateTime[end] - df.dateTime[1]).value ./ (1000 * 60 * 60 * 24)
println("N days: ", ndays)

pm_units = "μg⋅m⁻³"

# compute ticks for the figure

fmt_long = dateformat"m/yy"
Dates.format(df.dateTime[1], fmt_long)

# set up ticks to use minutes
mstart = df.dateTime[1] - Month(1)
mend = df.dateTime[end] + Month(1)

df.dt = [d.value / (1000*60) for d in df.dateTime .- mstart ]
df_bam.dt = [d.value / (1000*60) for d in df_bam.dateTime .- mstart ]

# Find the skips.
unique([ d.value ./ 1000 for d in (df.dateTime[2:end] .- df.dateTime[1:end-1])])

# compute skips
dt_val = 5*60
is_skip = [false, [v.value ./ (1000) != dt_val for v in df.dateTime[2:end] .- df.dateTime[1:end-1]]...]

# create list of data frames
idxs = findall(is_skip)
idxs = [1, idxs..., nrow(df)]
dfs = DataFrame[df[idxs[i]:idxs[i+1],:] for i in 1:length(idxs)-1];

for df in dfs
    println(nrow(df))
end


# Visualize the time series

pm1_0 = df.pm1_0;
pm2_5 = df.pm2_5;
pm10_0 = df.pm10_0;


# set up ticks to use
tick_values = [d.value / (1000*60) for d in collect(mstart:Month(1):mend) .- mstart]
minor_tick_values = [d.value / (1000*60) for d in collect(mstart:Day(7):mend) .- mstart]
tick_labels = [Dates.format(d, fmt_long) for d in mstart:Month(1):mend]
xticks = (tick_values, tick_labels)

# define colors for the PM time series

fig = Figure(size=fig_sizes.wide);
ax = Axis(
    fig[2,1],
    xlabel="time", xticks=xticks, xminorticks=minor_tick_values, xminorticksvisible=true, xminorgridvisible = false,
    ylabel="PM ($(pm_units))", yticks=0:20:200, yminorticksvisible=true, yminorgridvisible = false,
    );


ylims!(ax, 0, 100)
#ylims!(ax, 0, 20)
#ylims!(ax, 0, nothing)
xlims!(ax, df.dt[1], df.dt[end])

include_line = false

ls = []
for i in 1:length(dfs)
    b10  = band!(ax, dfs[i].dt, zeros(nrow(dfs[i])), dfs[i].pm10_0, color=(pm_colors[1]))
    b5   = band!(ax, dfs[i].dt, zeros(nrow(dfs[i])), dfs[i].pm5_0, color=(pm_colors[2]))
    b2_5 = band!(ax, dfs[i].dt, zeros(nrow(dfs[i])), dfs[i].pm2_5, color=(pm_colors[3]))
    b1   = band!(ax, dfs[i].dt, zeros(nrow(dfs[i])), dfs[i].pm1_0, color=(pm_colors[4]))
    b0_5 = band!(ax, dfs[i].dt, zeros(nrow(dfs[i])), dfs[i].pm0_5, color=(pm_colors[5]))
    b0_3 = band!(ax, dfs[i].dt, zeros(nrow(dfs[i])), dfs[i].pm0_3, color=(pm_colors[6]))
    b0_1 = band!(ax, dfs[i].dt, zeros(nrow(dfs[i])), dfs[i].pm0_1, color=(pm_colors[7]))

    if include_line
        #lines!(ax, dfs[i].dt, dfs[i].pm10_0, color=:black, linewidth=0.1)
        #lines!(ax, dfs[i].dt, dfs[i].pm5_0, color=:gray, linewidth=0.1)
        lines!(ax, dfs[i].dt, dfs[i].pm2_5, color=:darkgray, linewidth=0.1)
        #lines!(ax, dfs[i].dt, dfs[i].pm1_0, color=:lightgray, linewidth=0.1)
        #lines!(ax, dfs[i].dt, dfs[i].pm0_5, color=:gray, linewidth=0.1)
        #lines!(ax, dfs[i].dt, dfs[i].pm0_3, color=:gray, linewidth=0.1)
        #lines!(ax, dfs[i].dt, dfs[i].pm0_1, color=:gray, linewidth=0.1)
    end

    if i == 1
        push!(ls, b0_1)
        push!(ls, b0_3)
        push!(ls, b0_5)
        push!(ls, b1)
        push!(ls, b2_5)
        push!(ls, b5)
        push!(ls, b10)
    end
end

# add in floating legend to top of plot
fig[1,1] = Legend(
    fig,
    ls,
    [pm_labels["0.1"], pm_labels["0.3"], pm_labels["0.5"], pm_labels["1"], pm_labels["2.5"], pm_labels["5"], pm_labels["10"]],
    framevisible=false,
    orientation=:horizontal,
    labelsize=12,
    height=2, 
    halign=:center
) #, padding=(0,0,-18,0), labelsize=14, height=-5, halign=:center)

# xlims!(ax, minor_tick_values[10], minor_tick_values[12])
# xlims!(ax, tick_values[2], tick_values[3])

fig

save(joinpath(figpath, "1__pm-bands.png"), fig)
save(joinpath(figpath, "1__pm-bands.pdf"), fig)



# plot PM 1.0, 2.5, 10.0 time series

fig = Figure(size=fig_sizes.wide);
ax = Axis(
    fig[1,1],
    xlabel="time", xticks=xticks, xminorticks=minor_tick_values, xminorticksvisible=true, xminorgridvisible = false,
    #ylabel="PM 2.5 ($(pm_units))", yticks=0:20:200, yminorticksvisible=true, yminorgridvisible = false,
    ylabel=pm_labels_units["2.5"], yticks=0:20:200, yminorticksvisible=true, yminorgridvisible = false,
    );


ylims!(ax, 0, 100)
xlims!(ax, df.dt[1], df.dt[end])

for i in 1:length(dfs)
    l2_5 = lines!(ax, dfs[i].dt, dfs[i].pm2_5, color=mints_colors[3], linewidth=1, linecap=:round)
end

# xlims!(ax, minor_tick_values[10], minor_tick_values[12])
# xlims!(ax, tick_values[2], tick_values[3])

fig

save(joinpath(figpath, "1b_pm-2.5.png"), fig)
save(joinpath(figpath, "1b_pm-2.5.pdf"), fig)


# Visualize particle size distribution heatmap



size(PC)

y_edges = bin_edges["vec"]

crange = (-3, 7)
cmap = :thermal


fig = Figure(size=fig_sizes.wide);
ax = Axis(
    fig[1,1],
    xlabel="time", xticks=xticks, xminorticks=minor_tick_values, xminorticksvisible=true, xminorgridvisible = false,
    ylabel="Particle Diameter (μm)",
    yticks=0:1:10, ygridvisible=false, yminorgridvisible=false
    #yminorticks=y_edges, yminorticksvisible=true, yminorgridvisible=false
     );

maxes = []
mins = []

for i in 1:length(dfs)

    PC = Matrix(dfs[i][:, ["pc0_1", "pc0_3", "pc0_5", "pc1_0", "pc2_5", "pc5_0", "pc10_0"]])
    log10PC = log10.(PC)
    log10PC[isinf.(log10PC)] .= NaN

    #push!(maxes, maximum(log10PC[.!isnan.(log10PC)]))
    #push!(mins, minimum(log10PC[.!isnan.(log10PC)]))
    heatmap!(ax, dfs[i].dt, y_edges, log10PC, colorrange=crange, colormap=cmap)
end     

hlines!(ax, y_edges[2:end-1], color=:lightgray, linewidth=0.5, linestyle=:dash)


cb = Colorbar(fig[1,2], colorrange=crange, colormap=cmap, label=rich("log", subscript("10"), "(count⋅cm⁻³)"))

fig

save(joinpath(figpath, "1c_pc-heatmap.png"), fig)


# Evaluate assumed particle densities for each bin
PM  =  Matrix(dfs[1][:, ["pm0_1", "pm0_3", "pm0_5", "pm1_0", "pm2_5", "pm5_0", "pm10_0"]])
PC  =  Matrix(dfs[1][:, ["pc0_1", "pc0_3", "pc0_5", "pc1_0", "pc2_5", "pc5_0", "pc10_0"]])


size(PC)
size(PM)

coeffs = get_conversion_coeffs(dfs)


df_test = dfs[1][:, [:pc0_1, :pc0_3, :pc0_5, :pc1_0, :pc2_5, :pc5_0, :pc10_0]]

pm_from_pc!(df_test, coeffs)


msize = 7
al = 0.6

fig = Figure(size=fig_sizes.default);
ax = Axis(fig[1,1], xlabel="PM from IPS7100  ($(pm_units))", ylabel="PM from PC ($(pm_units))")

lines!(ax, [-500,500], [-500, 500], color=:darkgray, linewidth=2)
s10_0 = scatter!(ax, dfs[1].pm0_1, df_test.pm0_1, markersize=msize, alpha=al)
s5_0  = scatter!(ax, dfs[1].pm0_3, df_test.pm0_3, markersize=msize, alpha=al)
s2_5  = scatter!(ax, dfs[1].pm0_5, df_test.pm0_5, markersize=msize, alpha=al)
s1_0  = scatter!(ax, dfs[1].pm1_0, df_test.pm1_0, markersize=msize, alpha=al)
s0_5  = scatter!(ax, dfs[1].pm2_5, df_test.pm2_5, markersize=msize, alpha=al)
s0_3  = scatter!(ax, dfs[1].pm5_0, df_test.pm5_0, markersize=msize, alpha=al)
s0_1  = scatter!(ax, dfs[1].pm10_0, df_test.pm10_0, markersize=msize, alpha=al)

fig[1,2] = Legend(
    fig,
    [s0_1, s0_3, s0_5, s1_0, s2_5, s5_0, s10_0],
    [pm_labels["0.1"], pm_labels["0.3"], pm_labels["0.5"], pm_labels["1"], pm_labels["2.5"], pm_labels["5"], pm_labels["10"]],
    framevisible=false,
    orientation=:vertical,
    #labelsize=12,
    #height=2, 
    #halign=:center
) #, padding=(0,0,-18,0), labelsize=14, height=-5, halign=:center)

#
xlims!(ax, 0, 225)
ylims!(ax, 0, 225)

fig

save(joinpath(figpath, "2_pc-to-pm.png"), fig)
save(joinpath(figpath, "2_pc-to-pm.pdf"), fig)


# Now let's implement the correction

names(dfs[1])




coeffs = get_conversion_coeffs(dfs)

for df in dfs
    correct_pc!(df; kappa=0.62, eff=0.35)
    pm_from_pc_cor!(df, coeffs)
end


names(dfs[1])

dfs[1].pc1_0
dfs[1].pc1_0_cor

fig = Figure(size=fig_sizes.wide);

ax_top = Axis(
    fig[1,1], 
    xticks=xticks, xminorticks=minor_tick_values, xminorticksvisible=false, xminorgridvisible = true,
    xticksvisible=false, xticklabelsvisible=false, yminorgridvisible=false,
    ylabel="RH", yticks=0:25:100
)



ax = Axis(
    fig[2,1],
    xlabel="time", xticks=xticks, xminorticks=minor_tick_values, xminorticksvisible=true, xminorgridvisible = true,
    ylabel=pm_labels_units["2.5"], yticks=0:25:200, yminorticksvisible=true, yminorgridvisible = false,
);

linkxaxes!(ax, ax_top)
rowsize!(fig.layout, 1, Relative(0.25))

# l_bam = lines!(ax, df_bam.dt, df_bam.pm2_5BAM, color=:darkgray, linewidth=1)

ls = []
for i in 1:length(dfs)
    lines!(ax_top, dfs[i].dt, dfs[i].humidity, color=mints_colors[3], linewidth=1)
    l_orig = lines!(ax, dfs[i].dt, dfs[i].pm2_5, color=mints_colors[1], linewidth=1)
    l_cor = lines!(ax, dfs[i].dt, dfs[i].pm2_5_cor, color=mints_colors[2], linewidth=1)

    if i == 1
        push!(ls, l_orig)
        push!(ls, l_cor)
    end
end

# axislegend(ax, [ls..., l_bam], ["IPS7100", "IPS7100 HC", "BAM"], orientation=:horizontal, labelsize=10)
# axislegend(ax, ls, ["IPS7100", "IPS7100 HC"], orientation=:horizontal, labelsize=10)

# ylims!(ax, 0, 200)
# xlims!(ax, dfs[1].dt[1000], dfs[1].dt[2000])
# ylims!(ax, 0, 20)
xlims!(ax, dfs[1].dt[1], dfs[end].dt[end])

fig


