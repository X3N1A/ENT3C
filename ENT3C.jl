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

using Pkg
installs, resolve_env, config_file, julia_version = parse_args(ARGS)
#args = ["--config-file=config/config.json", "--install-deps=no"]
#installs, config_file = parse_args(args)

required_packages = ["DataFrames", "BenchmarkTools", "JSON", "Printf", "Plots", "ColorSchemes", "SuiteSparse", "HDF5", "NaNStatistics", "Statistics", "Combinatorics", "CSV"]
missing_packages = check_missing_packages(required_packages)

if isnothing(installs) && isnothing(resolve_env) && !isempty(missing_packages)
    println("Error! Missing julia packages. Please specify one of the following:")
    println("1. Use --install-deps=yes to install the packages globally.")
    println("   Or")
    println("2. Specify --resolve-env=yes with --julia-version=1.10.4 or 1.11.2 to install and resolve environment.")
    exit(1)
end

if installs == "yes" && length(missing_packages)!=1 && resolve_env != "yes"
    println("Installing missing dependencies globally.")
    println(resolve_env)
    install_packages(missing_packages)
end


if resolve_env == "yes"
        if isnothing(julia_version) || julia_version == ""
           throw(ArgumentError("--resolve-env requires julia version specification --julia-version=1.11.2 or --julia-version=1.10.4"))
        end
        println("Resolving ENT3C julia enviornment.")
        Pkg.activate("project_files/$(julia_version)")
        Pkg.add(["DataFrames", "BenchmarkTools", "JSON", "Printf", "Plots", "ColorSchemes", "SuiteSparse", "HDF5", "NaNStatistics", "Statistics", "Combinatorics", "CSV"])
        Pkg.instantiate()
end


using DataFrames, BenchmarkTools, JSON, Printf, Plots, ColorSchemes, SuiteSparse

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
weights_name = config["WEIGHTS_NAME"]

if OUT_PREFIX === nothing
    OUT_PREFIX = @sprintf("%dkb",Resolution/1e3)
end
OUT_DIR = config["OUT_DIR"]
OUT_DIR = string(OUT_DIR,"JULIA")
if !isdir(OUT_DIR)
    mkpath(OUT_DIR)
end

ChrNrs = string(config["ChrNr"])
if contains(ChrNrs, ",")
    ChrNrs = split(ChrNrs, ",")
else
    ChrNrs = [ChrNrs]
end
#############################################################################
# ENT3C table
##############################################################################
ENT3C_OUT = main(FILES,Resolutions,ChrNrs,SUB_M_SIZE_FIX,CHRSPLIT,PHI_MAX,phi,NormM,weights_name)
#f = [filter(row -> isnan(row.S)  , ENT3C_OUT).binNrStart[1], filter(row -> isnan(row.S)  , ENT3C_OUT).binNrEnd[1]]
#filter(row -> row.binNrStart==f[1] && row.binNrEnd==f[2], ENT3C_OUT)
#############################################################################
# similarity table
#############################################################################
SAMPLES::Array=unique(ENT3C_OUT.Name)
if length(SAMPLES)>1 
        Similarity = get_similarity_table(ENT3C_OUT,ChrNrs,Resolutions,Biological_replicates)
        CSV.write(@sprintf("%s/%s_ENT3C_similarity.csv",OUT_DIR,OUT_PREFIX), Similarity)
end
CSV.write(@sprintf("%s/%s_ENT3C_OUT.csv",OUT_DIR,OUT_PREFIX), ENT3C_OUT)
	
