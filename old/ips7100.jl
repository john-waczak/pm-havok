pm_units = "μg⋅m⁻³"

pm_colors = [
    colorant"#C3EE11",
    colorant"#2CD387",
    colorant"#2CCCD3",
    colorant"#2C78D3",
    colorant"#312CD3",
    colorant"#852ED1", 
    colorant"#D12ECC",
]



pm_labels_units = Dict(
    "0.1"   => rich("PM", subscript("0.1"), " ($(pm_units))"),
    "0.3"   => rich("PM", subscript("0.3"), " ($(pm_units))"),
    "0.5"   => rich("PM", subscript("0.5"), " ($(pm_units))"),
    "1"   => rich("PM", subscript("1"), " ($(pm_units))"),
    "2.5" => rich("PM", subscript("2.5"), " ($(pm_units))"),
    "5" => rich("PM", subscript("5"), " ($(pm_units))"),
    "10"  => rich("PM", subscript("10"), " ($(pm_units))"),

)


pm_labels = Dict(
    "0.1"   => rich("PM", subscript("0.1")),
    "0.3"   => rich("PM", subscript("0.3")),
    "0.5"   => rich("PM", subscript("0.5")),
    "1"   => rich("PM", subscript("1")),
    "2.5" => rich("PM", subscript("2.5")),
    "5" => rich("PM", subscript("5")),
    "10"  => rich("PM", subscript("10")),

)


bin_edges = Dict(
    "0.1" => (0, 0.1),
    "0.3" => (0.1, 0.3),
    "0.5" => (0.3, 0.5),
    "1" => (0.5, 1.0),
    "2.5" => (1.0, 2.5),
    "5" => (2.5, 5.0),
    "10" => (5.0, 10.0),
    "vec" => [0.0, 0.1, 0.3, 0.5, 1.0, 2.5, 5.0, 10.0],
)

bin_centers = (bin_edges["vec"][2:end] .- bin_edges["vec"][1:end-1])./2


# assume dfs is a list of dataframes where each has PC and PM components
function get_conversion_coeffs(dfs)
    As = []
    Bs = []
    Cs = []
    Ds = []
    Es = []
    Fs = [] 
    Gs = []

    for df in dfs
    
        A = mean([v for v in df.pm0_1 ./ df.pc0_1 if !isnan(v) && !isinf(v)])
        B = mean([v for v in (df.pm0_3 .- df.pm0_1) ./ df.pc0_3 if !isnan(v) && !isinf(v)])
        C = mean([v for v in (df.pm0_5 .- df.pm0_3) ./ df.pc0_5 if !isnan(v) && !isinf(v)])
        D = mean([v for v in (df.pm1_0 .- df.pm0_5) ./ df.pc1_0 if !isnan(v) && !isinf(v)])
        E = mean([v for v in (df.pm2_5 .- df.pm1_0) ./ df.pc2_5 if !isnan(v) && !isinf(v)])
        F = mean([v for v in (df.pm5_0 .- df.pm2_5) ./ df.pc5_0 if !isnan(v) && !isinf(v)])
        G = mean([v for v in (df.pm10_0 .- df.pm5_0) ./ df.pc10_0 if !isnan(v) && !isinf(v)])
    
        push!(As, A)
        push!(Bs, B)
        push!(Cs, C)
        push!(Ds, D)
        push!(Es, E)
        push!(Fs, F)
        push!(Gs, G)
    end

    A = mean([a for a in As if !isnan(a) && !isinf(a)])
    B = mean([b for b in Bs if !isnan(b) && !isinf(b)])
    C = mean([c for c in Cs if !isnan(c) && !isinf(c)])
    D = mean([d for d in Ds if !isnan(d) && !isinf(d)])
    E = mean([e for e in Es if !isnan(e) && !isinf(e)])
    F = mean([f for f in Fs if !isnan(f) && !isinf(f)])
    G = mean([g for g in Gs if !isnan(g) && !isinf(g)])

    coeffs = [A, B, C, D, E, F, G]
end


# assume df has no PM columns
function pm_from_pc!(df, coeffs)
    df.pm0_1  = coeffs[1] .* df.pc0_1
    df.pm0_3  = coeffs[2] .* df.pc0_3  .+ df.pm0_1
    df.pm0_5  = coeffs[3] .* df.pc0_5  .+ df.pm0_3
    df.pm1_0  = coeffs[4] .* df.pc1_0  .+ df.pm0_5
    df.pm2_5  = coeffs[5] .* df.pc2_5  .+ df.pm1_0
    df.pm5_0  = coeffs[6] .* df.pc5_0  .+ df.pm2_5
    df.pm10_0 = coeffs[7] .* df.pc10_0 .+ df.pm5_0
end



const DWET = [0.0, 0.1, 0.3, 0.5, 1.0, 2.5, 5.0, 10.0]

# RH    - relative humidity
# kappa - hygroscopic growth parameter
# eff   - Efflorescence point
function get_dry_bins(RH; kappa=0.62, eff=0.35, D_wet = DWET)
    rh = (RH == 100) ? 99.9 : RH
        
    if rh >= eff
        D_wet ./ ((1 + kappa * rh / (100.0 - rh))^(1/3))
    else
        return D_wet
    end
    
    # [(RH >= eff) ? D / (1 + kappa*RH/(100.0 - RH))^(1/3) : D for D in D_wet]
end


# b_low, b_high are corrected bin edges
function bin_fraction(b_low, b_high, B_low, B_high)
    @assert b_low <= b_high

    if b_low >= B_low && b_high <= B_high
        # new bin falls entirely within old bin
        return 1.0
    elseif b_low >= B_low && b_low < B_high && b_high > B_high        
        return (B_high - b_low)/(b_high - b_low)
    elseif b_low < B_low && b_high <= B_high && b_high > B_low
        return (b_high - B_low)/(b_high - b_low)
    elseif b_low < B_low && b_high > B_high
        # old bin falls entirely within new bin
        return (B_high-B_low)/(b_high-b_low) 
    elseif  b_high < B_low
        return 0.0
    elseif b_low > B_high
        return 0.0
    else
        println(b_low, "\t", b_low, "\t", B_low, "\t", B_high)
        return 0.0
    end
end


# given the counts in "correct" bin
# compute counts that go in original 
# bins
#
#
function bin_fractions(pc, b_low, b_high; D_wet = DWET)
    fracs = [
        bin_fraction(b_low, b_high, D_wet[1], D_wet[2]),
        bin_fraction(b_low, b_high, D_wet[2], D_wet[3]),
        bin_fraction(b_low, b_high, D_wet[3], D_wet[4]),
        bin_fraction(b_low, b_high, D_wet[4], D_wet[5]),
        bin_fraction(b_low, b_high, D_wet[5], D_wet[6]),
        bin_fraction(b_low, b_high, D_wet[6], D_wet[7]),
        bin_fraction(b_low, b_high, D_wet[7], D_wet[8]),
        ]

    return pc .* fracs
end




# assume df has PC and RH values
function correct_pc!(df; kappa=0.62, eff=0.35)
    # 1. add new columns to df
    df.pc0_1_cor = zeros(nrow(df))
    df.pc0_3_cor = zeros(nrow(df))
    df.pc0_5_cor = zeros(nrow(df))
    df.pc1_0_cor = zeros(nrow(df))
    df.pc2_5_cor = zeros(nrow(df))
    df.pc5_0_cor = zeros(nrow(df))
    df.pc10_0_cor = zeros(nrow(df))


    pc_0_1_out = zeros(7)
    pc_0_3_out = zeros(7)
    pc_0_5_out = zeros(7)
    pc_1_0_out = zeros(7)
    pc_2_5_out = zeros(7)
    pc_5_0_out = zeros(7)
    pc_10_0_out = zeros(7)

    # Loop over each row
    for row in eachrow(df)
        # zero out temp arrays
        pc_0_1_out .= 0.0
        pc_0_3_out .= 0.0
        pc_0_5_out .= 0.0
        pc_1_0_out .= 0.0
        pc_2_5_out .= 0.0
        pc_5_0_out .= 0.0
        pc_10_0_out .= 0.0
        
        
        # compute shifted bins
        D_dry = get_dry_bins(row.humidity, kappa=kappa, eff=eff)

        if D_dry[1] == D_dry[2]
            println("\n")
            println(D_dry)
            println(row)
            println("\n")
        end

        # redistribute values into new bins        
        pc_0_1_out  = bin_fractions(row.pc0_1,  D_dry[1], D_dry[2]) 
        pc_0_3_out  = bin_fractions(row.pc0_3,  D_dry[2], D_dry[3]) 
        pc_0_5_out  = bin_fractions(row.pc0_5,  D_dry[3], D_dry[4]) 
        pc_1_0_out  = bin_fractions(row.pc1_0,  D_dry[4], D_dry[5]) 
        pc_2_5_out  = bin_fractions(row.pc2_5,  D_dry[5], D_dry[6]) 
        pc_5_0_out  = bin_fractions(row.pc5_0,  D_dry[6], D_dry[7]) 
        pc_10_0_out = bin_fractions(row.pc10_0, D_dry[7], D_dry[8]) 


        row.pc0_1_cor  = pc_0_1_out[1] + pc_0_3_out[1] + pc_0_5_out[1] + pc_1_0_out[1] + pc_2_5_out[1] + pc_5_0_out[1] + pc_10_0_out[1]
        row.pc0_3_cor  = pc_0_1_out[2] + pc_0_3_out[2] + pc_0_5_out[2] + pc_1_0_out[2] + pc_2_5_out[2] + pc_5_0_out[2] + pc_10_0_out[2]
        row.pc0_5_cor  = pc_0_1_out[3] + pc_0_3_out[3] + pc_0_5_out[3] + pc_1_0_out[3] + pc_2_5_out[3] + pc_5_0_out[3] + pc_10_0_out[3]
        row.pc1_0_cor  = pc_0_1_out[4] + pc_0_3_out[4] + pc_0_5_out[4] + pc_1_0_out[4] + pc_2_5_out[4] + pc_5_0_out[4] + pc_10_0_out[4]
        row.pc2_5_cor  = pc_0_1_out[5] + pc_0_3_out[5] + pc_0_5_out[5] + pc_1_0_out[5] + pc_2_5_out[5] + pc_5_0_out[5] + pc_10_0_out[5]
        row.pc5_0_cor  = pc_0_1_out[6] + pc_0_3_out[6] + pc_0_5_out[6] + pc_1_0_out[6] + pc_2_5_out[6] + pc_5_0_out[6] + pc_10_0_out[6]
        row.pc10_0_cor = pc_0_1_out[7] + pc_0_3_out[7] + pc_0_5_out[7] + pc_1_0_out[7] + pc_2_5_out[7] + pc_5_0_out[7] + pc_10_0_out[7]                

    end
end



# assume df has pc_cor columns
function pm_from_pc_cor!(df, coeffs)
    df.pm0_1_cor  = coeffs[1] .* df.pc0_1_cor
    df.pm0_3_cor  = coeffs[2] .* df.pc0_3_cor  .+ df.pm0_1_cor
    df.pm0_5_cor  = coeffs[3] .* df.pc0_5_cor  .+ df.pm0_3_cor
    df.pm1_0_cor  = coeffs[4] .* df.pc1_0_cor  .+ df.pm0_5_cor
    df.pm2_5_cor  = coeffs[5] .* df.pc2_5_cor  .+ df.pm1_0_cor
    df.pm5_0_cor  = coeffs[6] .* df.pc5_0_cor  .+ df.pm2_5_cor
    df.pm10_0_cor = coeffs[7] .* df.pc10_0_cor .+ df.pm5_0_cor
end


