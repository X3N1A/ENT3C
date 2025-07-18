using CSV
using Combinatorics
using DataFrames
using HDF5
using LinearAlgebra
using NaNStatistics
using Plots
using Printf
using Statistics

struct INFO
    FN::String
    META::String # sample short name: e.g. G401_BR1
end
function parse_args(args)
    installs = nothing
    resolve_env = nothing
    config_file = nothing
    julia_version = nothing
    for arg in args
        if startswith(arg, "--config-file=")
            config_file = split(arg, "=")[2]
        end
        if startswith(arg, "--install-deps=")
            installs = split(arg, "=")[2]
        end
        if startswith(arg, "--resolve-env=")
            resolve_env = split(arg, "=")[2]
        end
        if startswith(arg, "--julia-version=")
            julia_version = split(arg, "=")[2]
        end
    end

    if isnothing(config_file) || config_file == ""
        throw(ArgumentError("--config-file is missing or empty"))
    end

    return installs, resolve_env, config_file, julia_version
end

function check_missing_packages(required_packages)
        installed_packages = [info.name for info in values(Pkg.dependencies())]
        missing_packages = [pkg for pkg in required_packages if !(pkg in installed_packages)]
        return missing_packages
end


function install_packages(packages)
        for pkg in packages
                try
                        Pkg.add(pkg)
                catch e
                        println("Failed to install package $pkg: $e")
                        exit(1)
                end
        end
end


function get_config(config_file)
    config=nothing
    if config_file === nothing
        println("Please specify a configuration file.")
    else
    println("Using config file: $config_file")
    try
        config = JSON.parsefile(config_file)
        println("Config data: ", config)
    catch e
        println("Error reading config file: ", e)
    end
    end
    return config
end

function check_input_option(varname, value, allowed_options)
    if !(isnothing(value) || value == "" || lowercase(value) in allowed_options)
        println(ArgumentError("ERROR: Illegal value for --$varname: '$value'. Allowed values are: $(collect(allowed_options))"))
        exit(1)
    end
end


function load_cooler(FN, ChrNr, Resolution, NormM, weights_name)
    # INPUT VARS
    # OUTPUT VARS
    preStr = ""
    if occursin("mcool", FN)
        preStr = "/resolutions/$Resolution"
    end
    # Load non-zero upper diagonal elements
    BINIDs = hcat(h5read(FN, "$preStr/pixels/bin1_id"), h5read(FN, "$preStr/pixels/bin2_id"))
    BINIDs = BINIDs .+ 1
    counts = h5read(FN, "$preStr/pixels/count")
    BINS = hcat(h5read(FN, "$preStr/bins/start"), h5read(FN, "$preStr/bins/end"))
    BINS = hcat(1:size(BINS, 1), BINS)
    if NormM == 1
        weights = h5read(FN, "$preStr/bins/$weights_name")
    else
        weights = fill(NaN, size(BINS, 1))
    end
  
    command = `h5dump -d $preStr/bins/chrom -H $FN`
    output = read(pipeline(command, `grep "^\s\+\""`, `sed -E 's/^\s*"([^"]+)"\s+([0-9]+);/\2 => "\1",/'`), String)
    output = """
    Dict(
    $(join(
        filter(line -> !isempty(line), split(output, '\n')),
        "\n"
    ))
    )
    """
    open("JULIA_functions/chrdict.jl", "w") do io
        write(io, output)
    end
    chrs = h5read(FN, "$preStr/bins/chrom") 
    chromosome_map = include("JULIA_functions/chrdict.jl")
    chrs = [chromosome_map[chr] for chr in chrs]

    BIN_TABLE = DataFrame(chrs=chrs, BINS_ALL=BINS[:, 1], START=BINS[:, 2], END=BINS[:, 3], weights=weights)
    BIN_TABLE = filter(row -> row.chrs == ChrNr || row.chrs == "chr$ChrNr" , BIN_TABLE)

    f = findall(
        (BINIDs[:, 1] .>= minimum(BIN_TABLE.BINS_ALL)) .&
        (BINIDs[:, 2] .>= minimum(BIN_TABLE.BINS_ALL)) .&
        (BINIDs[:, 1] .<= maximum(BIN_TABLE.BINS_ALL)) .&
        (BINIDs[:, 2] .<= maximum(BIN_TABLE.BINS_ALL)))

    BINIDs = BINIDs[f, :]
    counts = counts[f, :]
    BINIDs .= BINIDs .- BIN_TABLE.BINS_ALL[1] .+ 1

    CONTACT_MAP = fill(NaN, size(BIN_TABLE, 1), size(BIN_TABLE, 1))

    index = (BINIDs[:, 2] .- 1) .* size(CONTACT_MAP, 1) .+ BINIDs[:, 1]
    CONTACT_MAP[index] = counts

    if NormM == 1
        a = BIN_TABLE.weights
        b = BIN_TABLE.weights'
        CONTACT_MAP .= CONTACT_MAP .* (a*b)
    end
    CONTACT_MAP .= triu(CONTACT_MAP) .+ triu(CONTACT_MAP, 1)'
    INCLUDE = unique(findall(sum(isnan.(CONTACT_MAP), dims=2) .< size(CONTACT_MAP, 2)))
    EMPTY_BIN = fill(NaN, size(BIN_TABLE, 1))
    EMPTY_BIN[INCLUDE] .= 1

    BIN_TABLE[!, :CONTACT] = EMPTY_BIN
    BIN_TABLE[!, :binNr] = 1:size(BIN_TABLE, 1)

    return CONTACT_MAP, BIN_TABLE

#  BIN_TABLE[!,:binNr] .= 1:size(BIN_TABLE,1)
#  B = BIN_TABLE[.!isnan.(BIN_TABLE.CONTACT), :]
#  CONTACT_MAP = CONTACT_MAP[B.binNr,B.binNr]
#  CONTACT_MAP = replace!(CONTACT_MAP, NaN =>  minimum(vec(CONTACT_MAP[])))
#  gr()# pyplot(), gr(), or plotly()
#  heatmap(log.(CONTACT_MAP .+ 1e-10), aspect_ratio = 1, legend = false, color=:jet)

end


function vN_entropy(M::Matrix{Float64}, SUB_M_SIZE_FIX, CHRSPLIT, PHI_MAX, phi, BIN_TABLE, F)
    print(F, "\n")
    S = Vector{Float64}(undef, 0)
    N::Int = size(M, 1)

    if SUB_M_SIZE_FIX === nothing
        SUB_M_SIZE = round(N / CHRSPLIT)
    else
        SUB_M_SIZE = SUB_M_SIZE_FIX
    end

    PHI = Int(1 + floor((N - SUB_M_SIZE) ./ phi))
    if !isnothing(PHI_MAX)
        while PHI > PHI_MAX
            phi = phi + 1
            PHI = Int(1 + floor((N - SUB_M_SIZE) ./ phi))
        end
    end

    R1 = collect(1:phi:N)
    R1 = R1[1:PHI]
    R2 = R1 .+ SUB_M_SIZE .- 1
    R2 = R2[1:PHI]
    R = hcat(R1, R2)
    R = convert(Matrix{Int64}, R)

    replace!(M, 0.0 => NaN)
    M = log.(M)
    BIN_TABLE_NEW = Array{Int64}(undef, 0, 4)
    m = []
    for rr in 1:PHI
        m = M[R[rr, 1]:R[rr, 2], R[rr, 1]:R[rr, 2]]

        mask = isnan.(m) .| (m .== 0)
        mask =  sum(mask, dims=1)
 

        if any(vec(mask .>= (SUB_M_SIZE-1)))
            ENT = NaN 
        else
            replace!(m, NaN => minimum(filter(!isnan, m)))
            P = cor(m)

            SDs = std(m, dims=1) .< eps()
            P[SDs[1, :], :] .= 0
            P[:,SDs[1,:]] .= 0
            replace!(P, NaN => 0)
            rho = P ./ size(P, 1)

            #lam = eigen(rho);lam = lam.values
            lam = eigvals(rho)
            lam = lam[lam.>eps()]

            #gr()
            #heatmap(rho, aspect_ratio = 1, legend = false, color=:jet, axis = :square)
            ENT = -sum(real(lam .* log.(lam)))
        
        end
        S = vcat(S, ENT)

        BIN_TABLE_NEW = vcat(BIN_TABLE_NEW,
            [BIN_TABLE.binNr[R[rr, 1]], BIN_TABLE.binNr[R[rr, 2]],
                BIN_TABLE.START[R[rr, 1]], BIN_TABLE.END[R[rr, 2]]]')

    end

    return (S, SUB_M_SIZE, PHI, phi, BIN_TABLE_NEW)

end


function get_pairwise_combs(SAMPLES)

    comparisons = []
    SAMPLE_IDX = 1:length(SAMPLES)
    if length(SAMPLE_IDX) > 1
        comparisons = collect(combinations(SAMPLE_IDX, 2))

        comparisons = map(comp -> SAMPLES[comp], comparisons)
    else
        comparisons = SAMPLES
    end
    return comparisons

end

function get_similarity_table(ENT3C_OUT, ChrNrs, Resolutions, Biological_replicates, OUT_DIR,OUT_PREFIX)

    Similarity = DataFrame(Resolution=Vector{Int}[], ChrNr=Vector{Int}[], Sample1=Vector{String}[], Sample2=Vector{String}[], Q=Vector{Int}[])

    for Resolution in Resolutions
        PLT = []
        for ChrNr in ChrNrs
            SAMPLES::Array = unique(ENT3C_OUT.Name)
            comparisons = get_pairwise_combs(SAMPLES)
            pal = cgrad(:twelvebitrainbow)
            plt = plot()
            plotted = []
            c = 1
            for f in eachindex(comparisons)
                S1 = filter(row -> row.Name == comparisons[f][1] && row.ChrNr == ChrNr && row.Resolution == Resolution, ENT3C_OUT)
                S2 = filter(row -> row.Name == comparisons[f][2] && row.ChrNr == ChrNr && row.Resolution == Resolution, ENT3C_OUT)
                non_nan_idx = .!isnan.(S1.S).&.!isnan.(S2.S)
                Q = cor(S1.S[non_nan_idx], S2.S[non_nan_idx])

                Similarity = vcat(Similarity,
                    DataFrame(Resolution=Resolution, ChrNr=ChrNr, Sample1=comparisons[f][1],
                        Sample2=comparisons[f][2], Q=Q))

                if !(comparisons[f][1] in plotted)
                    plot!(plt, S1.S, linewidth=1, label=comparisons[f][1], linecolor=pal[c])
                    c = c + 1
                    xlims!(1, length(S1.S))
                    plotted = vcat(plotted, comparisons[f][1])
                end
                if !(comparisons[f][2] in plotted)
                    plot!(plt, S2.S, linewidth=1, label=comparisons[f][2], linecolor=pal[c])
                    c = c + 1
                    plotted = vcat(plotted, comparisons[f][2])
                end
            end
            plot!(plt, grid=true)
            if ChrNr != ChrNrs[end]
                plot!(plt, legend=false)
            else
                plot!(plt, legend=:outertopright)
            end
            #@df Similarity plot(1:size(Similarity,1), [:Sample1, :Sample2], group = :Replicate, colour = )
            if Biological_replicates == true
                cell_line(x) = match(r"^(.*?)_", x).captures[1] # r==regex, ^(.*?) = any chars, then underscore
                BR = filter(row -> cell_line(row.Sample1) == cell_line(row.Sample2) && row.ChrNr == ChrNr && row.Resolution == Resolution, Similarity)
                BR = mean(Float64.(BR.Q))
                print(BR)
                nonBR = filter(row -> cell_line(row.Sample1) != cell_line(row.Sample2) && row.ChrNr == ChrNr && row.Resolution == Resolution, Similarity)
                nonBR = mean(Float64.(nonBR.Q))
                title!(plt, @sprintf("Chr%s \$\\overline{Q}_{BR}=%4.2f\$ \$\\overline{Q}_{nonBR}=%4.2f\$", ChrNr, BR, nonBR), fontsize=12, interpreter=:latex)
            else
                title!(plt, @sprintf("Chr%s", ChrNr))
            end
            push!(PLT, plt)
        end
        PLT = plot(PLT..., size=(1801, 935))
        savefig(PLT, @sprintf("%s/%s_%d_ENT3C_OUT.svg", OUT_DIR, OUT_PREFIX, Resolution))
    end
    return (Similarity)
end


















