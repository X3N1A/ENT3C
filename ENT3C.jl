using DataFrames, BenchmarkTools, JSON, Printf

include("JULIA/ent3c_functions.jl")

config = JSON.parsefile("config/config.julia.json")

SUB_M_SIZE_FIX = config["SUB_M_SIZE_FIX"]
ChrNr = config["ChrNr"]
WN_MAX = config["WN_MAX"]
CHRSPLIT = config["CHRSPLIT"]
Resolution = config["Resolution"]
WS = config["WS"]
NormM = config["NormM"]
DATA_PATH = config["DATA_PATH"]
FILES = config["FILES"]
FILES = hcat(joinpath.(DATA_PATH, FILES[1:2:end-1]), FILES[2:2:end])
OUT_PREFIX = config["OUT_PREFIX"]
if OUT_PREFIX === nothing
    OUT_PREFIX = @sprintf("Chr%d_%dkb",ChrNr,Resolution/1e3)
end
OUT_DIR = config["OUT_DIR"]
if !isdir(OUT_DIR)
    mkpath(OUT_DIR)
end

ENT3C_OUT, Similarity, p = main(FILES,Resolution,ChrNr,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS)

savefig(p,@sprintf("%s/%s_ENT3C_OUT.png",OUT_DIR,OUT_PREFIX))
CSV.write(@sprintf("%s/%s_ENT3C_OUT.csv",OUT_DIR,OUT_PREFIX), ENT3C_OUT)
CSV.write(@sprintf("%s/%s_ENT3C_similarity.csv",OUT_DIR,OUT_PREFIX), Similarity)

