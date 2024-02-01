push!(LOAD_PATH,".")
using ENT3C, DataFrames, BenchmarkTools, Plots
#import .auxilary 
        # brings only the module name into scope; 
        # i need to access functions via e.g. auxilary.load_cooler; 
Nr_Matrices = 2

struct INFO
    FN::String
    META::String # sample short name: e.g. G401_BR1
end

FNs = [INFO("DATA_30e6/ENCSR079VIJ.BioRep1.mcool","G401_BR1"),
              INFO("DATA_30e6/ENCSR444WCZ.BioRep1.mcool","A549")]

ChrNr=19
Resolution::Int=40e3
SUB_M_SIZE_FIX=300
CHRSPLIT=7
WN_MAX::Int=1e3;
WS=1


OUT = DataFrame(Name = Vector{String}[], ChrNr = Vector{Int}[], Resolution = Vector{Int}[],
                n = Vector{Int}[], WN = Vector{Int}[], WS = Vector{Int}[],  binNrStart = Vector{Int}[], 
                binNrEnd = Vector{Int}[], START = Vector{Int}[], END = Vector{Int}[],  S = Vector{Float64}[])

#btime vN_entropy(M1,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS)
EXCLUDE =  Vector{Int64}(undef, 0)
for f in FNs
        
    FN = f.FN

    M, BIN_TABLE = load_cooler(FN, ChrNr, Resolution,0)

    EXCLUDE = vcat(EXCLUDE,BIN_TABLE.binNr[isnan.(BIN_TABLE.CONTACT).||isnan.(BIN_TABLE.weights)])
    EXCLUDE = unique(EXCLUDE)

    INCLUDE = collect(1:size(M,1))
    INCLUDE = setdiff(INCLUDE,EXCLUDE)
    M = M[INCLUDE,INCLUDE]

    BIN_TABLE = BIN_TABLE[INCLUDE,:]

    S, SUB_M_SIZE, WN, WS, BIN_TABLE_NEW = vN_entropy(M,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS,BIN_TABLE)

    N = length(INCLUDE)

    vcat(OUT,DataFrame(vec(repeat([f.META],N,1)),vec(repeat([ChrNr],N,1)),
         vec(repeat([Resolution],N,1)), vec(repeat([SUB_M_SIZE],N,1)),
         vec(repeat([WN],N,1)), vec(repeat([WS],N,1)),
         BIN_TABLE.START,BIN_TABLE.END,BIN_TABLE.START,BIN_TABLE.END,S))

    OUT = filter(row -> !(row.binNrStart in EXCLUDE)||!(row.binNrEnd in EXCLUDE), OUT)
end


















