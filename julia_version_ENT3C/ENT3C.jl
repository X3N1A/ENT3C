using Pkg,  DataFrames, BenchmarkTools, JSON, Printf, Plots, ColorSchemes, SuiteSparse

include("JULIA_functions/ent3c_functions.jl")


function main()
   
    installs, resolve_env, config_file, julia_version = parse_args(ARGS)
    config = get_config(config_file)

    #args = ["--config-file=config/config.json", "--install-deps=no","--julia-version=1.11.2"]
    #installs,  resolve_env, config_file, julia_version = parse_args(args)
    required_packages = ["DataFrames", "BenchmarkTools", "JSON", "Printf", "Plots", "ColorSchemes", "SuiteSparse", "HDF5", "NaNStatistics", "Statistics", "Combinatorics", "CSV"]
    missing_packages = check_missing_packages(required_packages)

    allowed_bool = Set(["yes", "true", "no", "false"])
    check_input_option("install-deps", installs, allowed_bool)
    check_input_option("resolve-env", resolve_env, allowed_bool)

    if !isnothing(julia_version) && julia_version != ""
        allowed_versions = Set(["1.10.4", "1.11.2"])
        if !(julia_version in allowed_versions)
            println(ArgumentError("ERROR: Illegal value for --julia-version: '$julia_version'. Allowed values are: $(collect(allowed_versions))"))
            exit(1)
        end
    end

    if isnothing(installs) && isnothing(resolve_env) && !isempty(missing_packages)
        println("ERROR: Missing julia packages. Please specify one of the following:")
        println("1. Use --install-deps=yes to install the packages globally.")
        println("   Or")
        println("2. Specify --resolve-env=yes with --julia-version=1.10.4 or 1.11.2 to install and resolve environment.")
        exit(1)
    end
    
    if (installs == "yes"||installs == "true") && length(missing_packages)!=1 && (resolve_env != "yes"||resolve_env != "true")
        println("Installing missing dependencies globally.")
        println(resolve_env)
        install_packages(missing_packages)
    end
    
    
    if (resolve_env == "yes"||resolve_env == "true")
            if isnothing(julia_version) || julia_version == ""
                println(ArgumentError("--resolve-env requires julia version specification --julia-version=1.11.2 or --julia-version=1.10.4"))
                exit(1)
            end
            println("Resolving ENT3C julia enviornment.")
            Pkg.activate("project_files/$(julia_version)")
            Pkg.add(["DataFrames", "BenchmarkTools", "JSON", "Printf", "Plots", "ColorSchemes", "SuiteSparse", "HDF5", "NaNStatistics", "Statistics", "Combinatorics", "CSV"])
            Pkg.instantiate()
    end    

    
	SUB_M_SIZE_FIX = config["SUB_M_SIZE_FIX"]
	PHI_MAX = config["PHI_MAX"]
	Resolutions = config["Resolution"]
	Resolutions = split(Resolutions, ",")
	Resolutions = [parse(Float64, res) |> Int for res in Resolutions]
	phi = config["phi"]
	NormM = config["NormM"]
	DATA_PATH = config["DATA_PATH"]
	FILES = config["FILES"]
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

    
    ENT3C_OUT = DataFrame(Name=Vector{String}[], ChrNr=Vector{Int}[], Resolution=Vector{Int}[],
        n=Vector{Int}[], PHI=Vector{Int}[], phi=Vector{Int}[], binNrStart=Vector{Int}[],
        binNrEnd=Vector{Int}[], START=Vector{Int}[], END=Vector{Int}[], S=Vector{Float64}[])


    for Resolution in Resolutions
        for ChrNr in ChrNrs
            FNs = []
            for f in 1:Int(size(FILES, 1))
                FNs = vcat(FNs, INFO(FILES[f, 1], FILES[f, 2]))
            end
 
            ##############################################################################
            # extract common empty bins of input matrices
            ##############################################################################
            EXCLUDE = 0
            for f in FNs
                FN = f.FN
                _, BIN_TABLE = load_cooler(FN, ChrNr, Resolution, NormM, weights_name)
                if NormM==0
                    EXCLUDE = vcat(EXCLUDE, BIN_TABLE.binNr[isnan.(BIN_TABLE.CONTACT)])
			    elseif NormM==1
				    EXCLUDE = vcat(EXCLUDE, BIN_TABLE.binNr[isnan.(BIN_TABLE.CONTACT)], BIN_TABLE.binNr[isnan.(BIN_TABLE.weights)])
			    end
            end
            EXCLUDE = unique(EXCLUDE)
            
            ##############################################################################
            #produce ENT3C_OUT data frame
            ##############################################################################
            for f in FNs
                
                FN = f.FN
                print(FN,"\n")
                M, BIN_TABLE = load_cooler(FN, ChrNr, Resolution, NormM, weights_name)

                INCLUDE = collect(1:size(M, 1))
                INCLUDE = setdiff(INCLUDE, EXCLUDE)
                M = M[INCLUDE, INCLUDE]

                BIN_TABLE = BIN_TABLE[INCLUDE, :]

                S, SUB_M_SIZE1, PHI_1, phi_1, BIN_TABLE_NEW = vN_entropy(M, SUB_M_SIZE_FIX, CHRSPLIT, PHI_MAX, phi, BIN_TABLE, FN)

                print(Resolution, " ", ChrNr, " ", FN, "\n")
                N = length(S)
                
                OUT1 = DataFrame(Name=fill(f.META, N), ChrNr=fill(ChrNr, N),
                    Resolution=fill(Resolution, N), n=fill(SUB_M_SIZE1, N),
                    PHI=fill(PHI_1, N), phi=fill(phi_1, N), binNrStart=BIN_TABLE_NEW[:, 1],
                    binNrEnd=BIN_TABLE_NEW[:, 2], START=BIN_TABLE_NEW[:, 3], END=BIN_TABLE_NEW[:, 4], S=S)

                ENT3C_OUT = vcat(ENT3C_OUT, OUT1)

                print(ENT3C_OUT,"\n")

            end
        end
    end


	#############################################################################
	# ENT3C table
	##############################################################################
	#f = [filter(row -> isnan(row.S)  , ENT3C_OUT).binNrStart[1], filter(row -> isnan(row.S)  , ENT3C_OUT).binNrEnd[1]]
	#filter(row -> row.binNrStart==f[1] && row.binNrEnd==f[2], ENT3C_OUT)
	#############################################################################
	# similarity table
	#############################################################################
	Biological_replicates = any(x -> contains(x, '_'), FILES[2, :])
	
	SAMPLES::Array=unique(ENT3C_OUT.Name)
	if length(SAMPLES)>1 
	        Similarity = get_similarity_table(ENT3C_OUT,ChrNrs,Resolutions,Biological_replicates, OUT_DIR, OUT_PREFIX)
	        CSV.write(@sprintf("%s/%s_ENT3C_similarity.csv",OUT_DIR,OUT_PREFIX), Similarity)
	end
	CSV.write(@sprintf("%s/%s_ENT3C_OUT.csv",OUT_DIR,OUT_PREFIX), ENT3C_OUT)
		
	return ENT3C_OUT, Similarity
end


ENT3C_OUT, Similarity = main()
