using HDF5, DataFrames, LinearAlgebra, Plots, NaNStatistics, Statistics, Combinatorics, Plots, CSV

struct INFO
        FN::String
        META::String # sample short name: e.g. G401_BR1
end

function main(FILES,Resolution,ChrNrs,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS,NormM)

        ENT3C_OUT = DataFrame(Name = Vector{String}[], ChrNr = Vector{Int}[], Resolution = Vector{Int}[],
	                        n = Vector{Int}[], WN = Vector{Int}[], WS = Vector{Int}[],  binNrStart = Vector{Int}[],
	                        binNrEnd = Vector{Int}[], START = Vector{Int}[], END = Vector{Int}[],  S = Vector{Float64}[])
        for ChrNr in ChrNrs
	
	        FNs=[]
	        for f in 1:Int(size(FILES,1))
	            FNs = vcat(FNs,INFO(FILES[f,1],FILES[f,2]))
	        end
	
	        ##############################################################################
	        # extract common empty bins of input matrices
	        ##############################################################################
	        EXCLUDE = 0
	        for f in FNs
	
	            FN = f.FN
	            M, BIN_TABLE = load_cooler(FN, ChrNr, Resolution,NormM)
	
	            EXCLUDE = vcat(EXCLUDE,BIN_TABLE.binNr[isnan.(BIN_TABLE.CONTACT)])
	        end
	        EXCLUDE = unique(EXCLUDE)
	        ##############################################################################
	        #produce ENT3C_OUT data frame
	        ##############################################################################
	        for f in FNs
	
	            FN = f.FN
	            M, BIN_TABLE = load_cooler(FN, ChrNr, Resolution,0)
	
	            INCLUDE = collect(1:size(M,1))
	            INCLUDE = setdiff(INCLUDE,EXCLUDE)
	            M = M[INCLUDE,INCLUDE]
	
	            BIN_TABLE = BIN_TABLE[INCLUDE,:]
	
	            S, SUB_M_SIZE1, WN1, WS1, BIN_TABLE_NEW = vN_entropy(M,SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS,BIN_TABLE)
	
	            N = length(S)
	
	            OUT1 = DataFrame(Name=fill(f.META,N),ChrNr=fill(ChrNr,N),
	                          Resolution = fill(Resolution,N), n = fill(SUB_M_SIZE1,N),
	                          WN = fill(WN1,N), WS = fill(WS1,N), binNrStart = BIN_TABLE_NEW[:,1],
	                          binNrEnd = BIN_TABLE_NEW[:,2],START = BIN_TABLE_NEW[:,3], END = BIN_TABLE_NEW[:,4],S =S)
	            ENT3C_OUT = vcat(ENT3C_OUT,OUT1)
	
	        end
	end
       return ENT3C_OUT
end


function load_cooler(FN, ChrNr, Resolution, NormM)
# INPUT VARS
# OUTPUT VARS
	preStr = ""
	if occursin("mcool", FN)
	        preStr = "/resolutions/$Resolution"
	end

        # Load non-zero upper diagonal elements
	BINIDs = hcat(h5read(FN, "$preStr/pixels/bin1_id"), h5read(FN, "$preStr/pixels/bin2_id"))
	counts = h5read(FN, "$preStr/pixels/count")
        
        BINS = hcat(h5read(FN, "$preStr/bins/start"), h5read(FN, "$preStr/bins/end"))
        BINS = hcat(1:size(BINS, 1), BINS)
        if NormM==1
                weights = h5read(FN, "$preStr/bins/weight")
            else
                 weights = fill(NaN, size(BINS,1))
        end

        chrs = h5read(FN, "$preStr/bins/chrom") .|> x -> x+1 .|> x -> "chr" .* string.(x)
        chrs[chrs .== "chr23"] .= "chrX"
        chrs[chrs .== "chr24"] .= "chrY"

	BIN_TABLE = DataFrame(chrs = chrs,BINS_ALL= BINS[:, 1],START = BINS[:, 2],END = BINS[:, 3],weights = weights)
        BIN_TABLE = filter(row ->  any(x -> "chr$x" == row.chrs, ChrNr), BIN_TABLE)
         
	f = findall(
                (BINIDs[:, 1] .>= minimum(BIN_TABLE.BINS_ALL)) .&
        	(BINIDs[:, 2] .>= minimum(BIN_TABLE.BINS_ALL)) .&
        	(BINIDs[:, 1] .<= maximum(BIN_TABLE.BINS_ALL)) .&
                (BINIDs[:, 2] .<= maximum(BIN_TABLE.BINS_ALL)))

	BINIDs = BINIDs[f, :]
	counts = counts[f, :]
	BINIDs .= BINIDs .- BINIDs[1, 1] .+ 1

	CONTACT_MAP = fill(NaN, size(BIN_TABLE, 1), size(BIN_TABLE, 1))

        index = (BINIDs[:, 2] .- 1) .* size(CONTACT_MAP, 1) .+ BINIDs[:, 1]
        CONTACT_MAP[index] = counts
        
        if  NormM==1
	        CONTACT_MAP .= CONTACT_MAP .* BIN_TABLE.weights
        end
	CONTACT_MAP .= triu(CONTACT_MAP) .+ triu(CONTACT_MAP, 1)'
       	INCLUDE = unique(findall(sum(isnan.(CONTACT_MAP), dims=2) .< size(CONTACT_MAP, 2)))
	EMPTY_BIN = fill(NaN, size(BIN_TABLE, 1))
	EMPTY_BIN[INCLUDE] .= 1

	BIN_TABLE[!, :CONTACT] = EMPTY_BIN
        BIN_TABLE[!, :binNr] = 1:size(BIN_TABLE,1) 

	return CONTACT_MAP, BIN_TABLE

        #BIN_TABLE[!,:binNr] .= 1:size(BIN_TABLE,1)
        #B = BIN_TABLE[.!isnan.(BIN_TABLE.CONTACT), :]
        #CONTACT_MAP = CONTACT_MAP[B.binNr,B.binNr]
        #CONTACT_MAP = CONTACT_MAP[600:800,600:800]
        #gr()# pyplot(), gr(), or plotly()
        #heatmap(log.(CONTACT_MAP), aspect_ratio = 1, legend = false, color=:jet, axis = :square)

end


function vN_entropy(M::Matrix{Float64},SUB_M_SIZE_FIX,CHRSPLIT,WN_MAX,WS,BIN_TABLE)
        S=Vector{Float64}(undef, 0)
        N::Int=size(M,1)
	if SUB_M_SIZE_FIX==nothing
	    SUB_M_SIZE=round(N/CHRSPLIT)
        else
            SUB_M_SIZE = SUB_M_SIZE_FIX
        end

        WN=Int(1+floor((N-SUB_M_SIZE)./WS))
        if WN_MAX!=nothing
            while WN>WN_MAX
               WS=WS+1
               WN=Int(1+floor((N-SUB_M_SIZE)./WS))
            end
        end


        R1=collect(1:WS:N)
        R1=R1[1:WN]
        R2=R1.+SUB_M_SIZE.-1
        R2=R2[1:WN]
        R=hcat(R1,R2)
        R=convert(Matrix{Int64}, R)

        M=log.(M)
        BIN_TABLE_NEW = Array{Int64}(undef, 0,4)
	for rr in 1:WN
            m=M[R[rr,1]:R[rr,2],R[rr,1]:R[rr,2]]
            replace!(m,NaN=>minimum(filter(!isnan,m)))
            P = cor(m);

            zero_variance_indices = all(m .== 0, dims=2) .& all(m .== 0, dims=1)
            P[zero_variance_indices] .= 0
            replace!(P,NaN=>0)

            rho=P./size(P,1);
		
            #lam = eigen(rho);lam = lam.values
            lam = eigvals(rho)
	    lam = lam[lam.>0];
	
            #gr()
            #heatmap(rho, aspect_ratio = 1, legend = false, color=:jet, axis = :square)


	    ENT = -sum(real(lam.*log.(lam)));
	
	    S = vcat(S,ENT)

            BIN_TABLE_NEW = vcat(BIN_TABLE_NEW,
                                 [BIN_TABLE.binNr[R[rr,1]],BIN_TABLE.binNr[R[rr,2]],
                                  BIN_TABLE.START[R[rr,1]],BIN_TABLE.END[R[rr,2]]]')
	end

	return(S,SUB_M_SIZE,WN,WS,BIN_TABLE_NEW)
	
end


function get_pairwise_combs(SAMPLES)


        SAMPLE_IDX = 1:length(SAMPLES)
        
        comparisons = collect(combinations(SAMPLE_IDX,2))

        comparisons = map(comp -> SAMPLES[comp], comparisons)

        return comparisons

end

function get_similarity_table(ENT3C_OUT,Biological_replicates)

        Similarity = DataFrame(ChrNr= Vector{Int}[], Sample1 = Vector{String}[], Sample2 = Vector{String}[],Q = Vector{Int}[])
        PLT=[]
        for ChrNr in ChrNrs
	        SAMPLES::Array=unique(ENT3C_OUT.Name)
	        comparisons =  get_pairwise_combs(SAMPLES)
                pal = cgrad(:twelvebitrainbow)
	        plt = plot()
	        plotted = []
	        c=1
	        for f in 1:size(comparisons,1)

	                S1 = filter(row -> row.Name==comparisons[f][1]&&row.ChrNr==ChrNr,ENT3C_OUT)
	                S2 = filter(row -> row.Name==comparisons[f][2]&&row.ChrNr==ChrNr,ENT3C_OUT)
	                Q = cor(S1.S,S2.S)
	
	                Similarity = vcat(Similarity,
	                DataFrame(ChrNr = ChrNr, Sample1 = comparisons[f][1],
	                Sample2 = comparisons[f][2], Q = Q))
	
	                if !(comparisons[f][1] in plotted)
	                        plot!(plt,S1.S, linewidth=1,label=comparisons[f][1], linecolor=pal[c]);c=c+1
	                        xlims!(1,length(S1.S))
	                        plotted = vcat(plotted,comparisons[f][1])
	                end
	                if !(comparisons[f][2] in plotted)
	                        plot!(plt,S2.S, linewidth=1,label=comparisons[f][2], linecolor=pal[c]);c=c+1
                                plotted = vcat(plotted,comparisons[f][2])
	                end
                end
	        plot!(plt, grid=true)
                if ChrNr!=ChrNrs[end]
	                plot!(plt, legend=false)
                else
                        plot!(plt, legend=:outertopright)
                end
	        #@df Similarity plot(1:size(Similarity,1), [:Sample1, :Sample2], group = :Replicate, colour = )
                if Biological_replicates==true
        	        cell_line(x) = match(r"^(.*?)_?", x).captures[1] # r==regex, ^(.*?) = any chars, then underscore
	                BR = filter(row -> cell_line(row.Sample1)==cell_line(row.Sample2) && row.ChrNr==ChrNr,Similarity)
	                BR = mean(BR.Q)
	                nonBR = filter(row -> cell_line(row.Sample1)!=cell_line(row.Sample2) && row.ChrNr==ChrNr,Similarity)
	                nonBR = mean(nonBR.Q)
	                title!(plt,@sprintf("Chr%d \$\\overline{Q}_{BR}=%4.2f\$ \$\\overline{Q}_{nonBR}=%4.2f\$", ChrNr, BR, nonBR), fontsize=12, interpreter=:latex)
                else
                        title!(plt,@sprintf("Chr%d",ChrNr))
                end
	        push!(PLT, plt)
        end
        PLT=plot(PLT...,size=(1801,935))

        return(Similarity,PLT)
end


















