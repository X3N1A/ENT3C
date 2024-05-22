using DataFrames, BenchmarkTools, JSON, Printf, Plots, ColorSchemes

include("JULIA_functions/ent3c_functions.jl")

config = JSON.parsefile("config/config.pooledBRs.json")

SUB_M_SIZE_FIX = config["SUB_M_SIZE_FIX"]
WN_MAX = config["WN_MAX"]
Resolutions = config["Resolution"]
Resolutions = split(Resolutions, ",")
Resolutions = [parse(Float64, res) |> Int for res in Resolutions]
WS = config["WS"]
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

if contains(config["ChrNr"],":")
        ChrNrs_1, ChrNrs_2 = parse.(Int, split(config["ChrNr"], ':'))
        ChrNrs = ChrNrs_1:ChrNrs_2
    else
        ChrNrs = parse.(Int, config["ChrNr"])
end

#############################################################################
# ENT3C table
##############################################################################
ENT3C_OUT = main(FILES,Resolutions,ChrNrs,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS,NormM)
#ENT3C_OUT =  CSV.File(@sprintf("%s/%s_ENT3C_OUT.csv",OUT_DIR,OUT_PREFIX)) |> DataFrame

#############################################################################
# similarity table
##############################################################################
Similarity = get_similarity_table(ENT3C_OUT,ChrNrs,Resolutions,Biological_replicates)

CSV.write(@sprintf("%s/%s_ENT3C_OUT.csv",OUT_DIR,OUT_PREFIX), ENT3C_OUT)
CSV.write(@sprintf("%s/%s_ENT3C_similarity.csv",OUT_DIR,OUT_PREFIX), Similarity)
	
