using DataFrames, BenchmarkTools, JSON, Printf, Plots, ColorSchemes

include("JULIA_functions/ent3c_functions.jl")

config = JSON.parsefile("config/config.multiResolution.pooledBRs.json")

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
    ChrNrs = [parse(Float64, c) |> Int for c in ChrNrs]
elseif contains(ChrNrs, ":")
    ChrNrs = split(ChrNrs, ":")
    ChrNrs = [parse(Float64, c) |> Int for c in ChrNrs]
    ChrNrs = collect(ChrNrs[1]:ChrNrs[end])
elseif isinteger(parse(Int8, ChrNrs))
    ChrNrs = parse(Int8, ChrNrs)
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
	
