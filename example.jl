#!/usr/local/JULIA/julia-1.10.0/bin/julia
push!(LOAD_PATH,".")
using DataFrames, BenchmarkTools, Plots, Statistics, CSV

include("ent3c_functions.jl")

struct INFO
    FN::String
    META::String # sample short name: e.g. G401_BR1
end

FNs = [INFO("DATA_30e6/ENCSR079VIJ.BioRep1.mcool","G401_BR1"),
       INFO("DATA_30e6/ENCSR079VIJ.BioRep2.mcool","G401_BR2"),
       INFO("DATA_30e6/ENCSR444WCZ.BioRep1.mcool","A549_BR1"),
       INFO("DATA_30e6/ENCSR444WCZ.BioRep2.mcool","A549_BR2")]
	
function ENT3C(FNs,Resolution,ChrNr,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS)

	#btime vN_entropy(M1,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS)
	##############################################################################
	#need to filter matrices first otherwise the windows wont coincide in final DF
	##############################################################################
	EXCLUDE = 0
	for f in FNs
	        
	    FN = f.FN
	    M, BIN_TABLE = load_cooler(FN, ChrNr, Resolution,0)
	    
	    EXCLUDE = vcat(EXCLUDE,BIN_TABLE.binNr[isnan.(BIN_TABLE.CONTACT).||isnan.(BIN_TABLE.weights)])
	
	end
	
	EXCLUDE = unique(EXCLUDE)
	
	##############################################################################
	#get ENT3C data frame
	##############################################################################
	
	ENT3C_OUT = DataFrame(Name = Vector{String}[], ChrNr = Vector{Int}[], Resolution = Vector{Int}[],
	                n = Vector{Int}[], WN = Vector{Int}[], WS = Vector{Int}[],  binNrStart = Vector{Int}[], 
	                binNrEnd = Vector{Int}[], START = Vector{Int}[], END = Vector{Int}[],  S = Vector{Float64}[])
	
	for f in FNs
	
	    FN = f.FN
	    M, BIN_TABLE = load_cooler(FN, ChrNr, Resolution,0)
	
	    INCLUDE = collect(1:size(M,1))
	    INCLUDE = setdiff(INCLUDE,EXCLUDE)
	    M = M[INCLUDE,INCLUDE]
	
	    BIN_TABLE = BIN_TABLE[INCLUDE,:]
	
	    S, SUB_M_SIZE1, WN1, WS1, BIN_TABLE_NEW = vN_entropy(M,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS,BIN_TABLE)
	
	    N = length(S)
	    print(FN)
	
	    OUT1 = DataFrame(Name=fill(f.META,N),ChrNr=fill(ChrNr,N),
	                  Resolution = fill(Resolution,N), n = fill(SUB_M_SIZE1,N),
	                  WN = fill(WN1,N), WS = fill(WS1,N), binNrStart = BIN_TABLE_NEW[:,1],
	                  binNrEnd = BIN_TABLE_NEW[:,2],START = BIN_TABLE_NEW[:,3], END = BIN_TABLE_NEW[:,4],S =S)
	    ENT3C_OUT = vcat(ENT3C_OUT,OUT1)
	
	end
	
	
	#############################################################################
	#get ENT3C data frame
	##############################################################################
	
	SAMPLES::Array=unique(ENT3C_OUT.Name)
	comparisons =  get_pairwise_combs(SAMPLES)
	
	Similarity = DataFrame(Sample1 = Vector{String}[], Sample2 = Vector{String}[],Q = Vector{Int}[])
	plotted = []
	p = plot()
	for f in 1:size(comparisons,1)
	        
	        S1 = filter(row -> row.Name==comparisons[f][1],ENT3C_OUT)
	        S2 = filter(row -> row.Name==comparisons[f][2],ENT3C_OUT)
	        Q = cor(S1.S,S2.S)
	    
	        Similarity = vcat(Similarity,
	                          DataFrame(Sample1 = comparisons[f][1],
	                                    Sample2 = comparisons[f][2], Q = Q))
	
	        if !(comparisons[f][1] in plotted)
	                plot!(p,S1.S, label=comparisons[f][1], linewidth=3)
	                plotted = vcat(plotted,comparisons[f][1])
	        end
	        if !(comparisons[f][2] in plotted)
	                plot!(p,S2.S, label=comparisons[f][2], linewidth=3)
	                plotted = vcat(plotted,comparisons[f][2])
	        end
	end
	
	# display(p)
	
        return(ENT3C_OUT,Similarity,p)	
end


Resolution::Int=40e3
ChrNr::Int=22
SUB_M_SIZE_FIX::Int=0
CHRSPLIT::Int=7
WN_MAX::Int=1e3
WS::Int=1
	
ENT3C_OUT, Similarity, p = ENT3C(FNs,Resolution,ChrNr,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS)

savefig(p,"ENT3C_OUT.png")
CSV.write("ENT3C_OUT.csv", ENT3C_OUT)
CSV.write("ENT3C_similarity.csv", Similarity)










