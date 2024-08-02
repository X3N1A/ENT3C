function parse_args(args)
    installs = nothing
    config_file = nothing
    for arg in args
        if startswith(arg, "--config-file=")
            config_file = split(arg, "=")[2]
        end
        if startswith(arg, "--install-deps=")
            installs = split(arg, "=")[2]
        end
    end
    return installs, config_file
end

installs, config_file = parse_args(ARGS)
#installs, config_file = parse_args(args)

if cmp(installs,"yes") == 0
        println("Installing missing dependencies and resolving enviornment.")
        using Pkg

        Pkg.activate(".")
        Pkg.instantiate()
end

using DataFrames, BenchmarkTools, JSON, Printf, Plots, ColorSchemes

include("JULIA_functions/ent3c_functions.jl")
######################################################
# user input config file
######################################################
#config = JSON.parsefile("config/config.test.json")

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
config = get_config(config_file)
#########################################################

SUB_M_SIZE_FIX = config["SUB_M_SIZE_FIX"]
PHI_MAX = config["PHI_MAX"]
Resolutions = config["Resolution"]
Resolutions = split(Resolutions, ",")
Resolutions = [parse(Float64, res) |> Int for res in Resolutions]
phi = config["phi"]
NormM = config["NormM"]
DATA_PATH = config["DATA_PATH"]
FILES = config["FILES"]
Biological_replicates = any(x -> contains(x, '_'), FILES[2, :])
FILES = hcat(joinpath.(DATA_PATH, FILES[1:2:end-1]), FILES[2:2:end])
OUT_PREFIX = config["OUT_PREFIX"]
CHRSPLIT = config["CHRSPLIT"]

if OUT_PREFIX === nothing
    OUT_PREFIX = @sprintf("%dkb",Resolution/1e3)
end
OUT_DIR = config["OUT_DIR"]
OUT_DIR = string(OUT_DIR,"JULIA")
if !isdir(OUT_DIR)
    mkpath(OUT_DIR)
end

ChrNrs = config["ChrNr"]
if contains(ChrNrs, ",")
    ChrNrs = split(ChrNrs, ",")
end
#############################################################################
# ENT3C table
##############################################################################
ENT3C_OUT = main(FILES,Resolutions,ChrNrs,SUB_M_SIZE_FIX,CHRSPLIT,PHI_MAX,phi,NormM)

#############################################################################
# similarity table
#############################################################################
SAMPLES::Array=unique(ENT3C_OUT.Name)
if length(SAMPLES)>1 
        Similarity = get_similarity_table(ENT3C_OUT,ChrNrs,Resolutions,Biological_replicates)
        CSV.write(@sprintf("%s/%s_ENT3C_similarity.csv",OUT_DIR,OUT_PREFIX), Similarity)
end
CSV.write(@sprintf("%s/%s_ENT3C_OUT.csv",OUT_DIR,OUT_PREFIX), ENT3C_OUT)
	
