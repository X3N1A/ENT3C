using HDF5, DataFrames, LinearAlgebra, Plots

function load_cooler(FN, ChrNr, Resolution)
# INPUT VARS
        FN="/HDD3/DocumentsHDD3/07_2022-/Mapping/ENCSR444WCZ_A549/PreProcessed/BioRep2/ENCSR444WCZ.BioRep2.30000000.40kb.cool"
# OUTPUT VARS
   # microC ... balanced matrix at desired Resolution
   # BIN_TABLE ... DataFrame containing bin information        

# Check if the file name contains 'mcool'
	preStr = ""
	if occursin("mcool", FN)
	        preStr = "/resolutions/$Resolution"
	end

# Load non-zero upper diagonal elements
	BINIDs = hcat(h5read(FN, "$preStr/pixels/bin1_id"), h5read(FN, "$preStr/pixels/bin2_id"))
	counts = h5read(FN, "$preStr/pixels/count")
        
        BINS = hcat(h5read(FN, "$preStr/bins/start"), h5read(FN, "$preStr/bins/end"))
        BINS = hcat(1:size(BINS, 1), BINS)
        weights = h5read(FN, "$preStr/bins/weight")

        chrs = h5read(FN, "$preStr/bins/chrom") .|> x -> x+1 .|> x -> "chr" .* string.(x)
        chrs[chrs .== "chr23"] .= "chrX"
        chrs[chrs .== "chr24"] .= "chrY"

	BIN_TABLE = DataFrame(chrs = chrs,binNrALL = BINS[:, 1],
	        START = BINS[:, 2],END = BINS[:, 3],weights = weights)
        
        
        BIN_TABLE = filter(row ->  any(x -> "chr$x" == row.chrs, ChrNr), BIN_TABLE)
        
	f = findall(
        	(BINIDs[:, 1] .>= minimum(BIN_TABLE.binNrALL)) .&
        	(BINIDs[:, 2] .>= minimum(BIN_TABLE.binNrALL)) .&
        	(BINIDs[:, 1] .<= maximum(BIN_TABLE.binNrALL)) .&
                (BINIDs[:, 2] .<= maximum(BIN_TABLE.binNrALL)))

	BINIDs = BINIDs[f, :]
	counts = counts[f, :]
	BINIDs .= BINIDs .- BINIDs[1, 1] .+ 1

	CONTACT_MAP = fill(NaN, size(BIN_TABLE, 1), size(BIN_TABLE, 1))

        index = (BINIDs[:, 2] .- 1) .* size(CONTACT_MAP, 1) .+ BINIDs[:, 1]
        CONTACT_MAP[index] = counts

	CONTACT_MAP .= CONTACT_MAP .* BIN_TABLE.weights
	CONTACT_MAP .= triu(CONTACT_MAP) .+ triu(CONTACT_MAP, 1)'
       	INCLUDE = unique(findall(sum(isnan.(CONTACT_MAP), dims=2) .< size(CONTACT_MAP, 2)))
	EMPTY_BIN = fill(NaN, size(BIN_TABLE, 1))
	EMPTY_BIN[INCLUDE] .= 1
	EMPTY_BIN[isnan.(BIN_TABLE.weights)] .= NaN

	BIN_TABLE[!, :CONTACT] = EMPTY_BIN

	return CONTACT_MAP, BIN_TABLE

        #BIN_TABLE[!,:binNrCHRS] .= 1:size(BIN_TABLE,1)
        #B = BIN_TABLE[.!isnan.(BIN_TABLE.CONTACT), :]
        #CONTACT_MAP = CONTACT_MAP[B.binNrCHRS,B.binNrCHRS]
        #CONTACT_MAP = CONTACT_MAP[600:800,600:800]
        #gr()# pyplot(), gr(), or plotly()
        #heatmap(log.(CONTACT_MAP), aspect_ratio = 1, legend = false, color=:jet, axis = :square)


end
