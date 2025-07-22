<figure>
    <img src="Figures/LogoSmall.png" width="200" >
</figure>

# 

ENT3C is a method for qunatifying the similarity of micro-C/Hi-C derived chromosomal contact matrices. It is based on the von Neumann entropy<sup>1</sup> and recent work on entropy quantification of Pearson correlation matrices<sup>2</sup>.
For a contact matrix, ENT3C records the change in local pattern *complexity* of smaller Pearson-transformed submatrices along a matrix diagonal to generate a characteristic signal. Similarity is defined as the Pearson correlation between the respective entropy signals of two contact matrices.

https://doi.org/10.1093/nargab/lqae076

## Summary of ENT3C approach
1. Loads cooler files and looks for shared empty bins.
2. ENT3C will first take the logarithm of an input matrix $\mathbf{M}$
2. Next, smaller submatrices $\mathbf{a}$ of dimension $n\times n$ are extracted along the diagonal of an input contact matrix $\mathbf{M}$
4. $nan$ values in $\mathbf{a}$ are set to the minimum value in $\mathbf{a}$.
5. $\mathbf{a}$ is transformed into a Pearson correlation matrix $\mathbf{P}$.
6. $\mathbf{P}$ is transformed into $\boldsymbol{\rho}=\mathbf{P}/n$ to fulfill the conditions for computing the von Neumann entropy.
7. The von Neumann entropy of $\boldsymbol{\rho}$ is computed as

   $S(\boldsymbol{\rho})=\sum_j \lambda_j \log \lambda_j$

   where $\lambda_j$ is the $j$ th eigenvalue of $\boldsymbol{\rho}$
8. This is repeated for subsequent submatrices along the diagonal of the input matrix and stored in the *entropy signal* $\mathbf{S}\_{M}$.
9. Similarity $Q$ is defined as the Pearson correlation $r$ between the entropy signals of two matrices: $Q(\mathbf{M}\_1,\mathbf{M}\_2) = r(\mathbf{S}\_{\mathbf{M}\_1},\mathbf{S}\_{\mathbf{M}\_2})$.

<figure>
    <img src="Figures/ENT3C_explain_2cells.png" width="400" 
         alt="explaination of ENT3C">
</figure>

Exemplary epiction of ENT3C derivation of the entropy signal $\mathbf{S}$ of two contact matrices $\mathbf{M}\_1$ and $\mathbf{M}\_2$. ENT3C's was run with  submatrix dimension $n=300$, window shift $\varphi=10$, and maximum number of data points in $\boldsymbol{S}$, $\Phi\_{\max}=\infty$, resulting in $\Phi=147$ submatrices. For subsequent scaled Pearson-transformed submatrices, $\boldsymbol{\rho}\_i$, along the diagonal of $\log{\boldsymbol{M}}$, ENT3C computes the von Neumann entropies $S(\boldsymbol{\rho}\_1), S(\boldsymbol{\rho}\_2), \ldots, S(\boldsymbol{\rho}\_{\Phi})$. The resulting signal $\mathbf{S} = \langle S(\boldsymbol{\rho}\_{1}), S(\boldsymbol{\rho}\_{2}), \ldots, S(\boldsymbol{\rho}\_{\Phi}) \rangle$ is shown in blue under the matrix. The first two ($\boldsymbol{\rho}\_{1-2}$), middle ($\boldsymbol{\rho}\_{73}$), and last two submatrices ($\boldsymbol{\rho}\_{146-147}$) are shown.


# Requirements and Installation

### Python (>=3.11), Julia or MATLAB. 

### Python:

*  https://pypi.org/project/ENT3C/

1) generate and activate python environment 
	
	```
	python3.12 -m venv .ent3c\_venv

	source .ent3c\_venv/bin/activate
	```

2) install ENT3C and requirements via ```pyproject.toml```: 

	```
	pip install .
	```

  
### Julia:

* packages: DataFrames, BenchmarkTools, JSON, Printf, Plots, ColorSchemes, SuiteSparse, HDF5, NaNStatistics, Statistics, Combinatorics, CSV
* For the Julia implementation, ubuntu's hdf5-tools is also required 
* Initial julia set-up
	
	* option for automatic global installation ```--install-deps=yes```. (Works with any julia version)
	* predefined julia enviornments for julia versions 1.10.4 or 1.11.2 are defined in ```project_files/<v.v.v>/Manifest.toml``` and ```project_files/<v.v.v>/Project.toml```
	* option to load enviornments with ```--resolve-env=yes``` and ```--julia-version=<v.v.v>```



# Parameters and configuration files of ENT3C

* all ENT3C parameters are defined in .json files ```config/config.json```

* The main ENT3C parameter affecting the final entropy signal $S$ is the dimension of the submatrices ```SUB_M_SIZE_FIX```. 

	* ```SUB_M_SIZE_FIX``` can be either be fixed by or alternatively, one can specify ```CHRSPLIT```; in this case ```SUB_M_SIZE_FIX``` will be computed internally to fit the number of desired times the contact matrix is to be paritioned into. 

	  ```PHI=1+floor((N-SUB_M_SIZE)./phi)```

	  where ```N``` is the size of the input contact matrix, ```phi``` is the window shift, ```PHI``` is the number of evaluated submatrices (consequently the number of data points in $S$).

* All implementations (```ENT3C.py```, ```ENT3C.jl``` and ```ENT3C.m```) use a configuration file in JSON format. 
	* example can be found in <config/config.json>


1) ```"DATA_PATH": "DATA"``` $\dots$ input data path. 

2) input files in format: ```[<COOL_FILENAME>, <SHORT_NAME>]```
``` 
"FILES": [
	"ENCSR079VIJ.BioRep1.40kb.cool",
 
	"G401_BR1",
 
	"ENCSR079VIJ.BioRep2.40kb.cool",
 
	"G401_BR2"]
``` 
- ENT3C also takes ```mcool``` files as input. Please refer to biological replicates as "_BR%d" in the <SHORT_NAME>.

&#9888; if comparing biological replicate samples, please ensure they are indicated as <_BR\#> in the config file &#9888;

4) ```"`OUT_DIR": "OUTPUT/"``` $\dots$ output directory. ```OUT_DIR``` will be concatenated with ```OUTPUT/JULIA/``` or ```OUTPUT/MATLAB/```.

5) ```"OUT_PREFIX": "40kb"``` $\dots$ prefix for output files.

6) ```"Resolution": "40e3,100e3"``` $\dots$ resolutions to be evaluated. 

7) ```"ChrNr": "15,16,17,18,19,20,21,22,X"``` $\dots$ chromosome numbers to be evaluated.

8) ```"NormM": 0``` $\dots$ input contact matrices can be balanced. If ```NormM: 1```, balancing weights in cooler are applied. If set to 1, ENT3C expects weights to be in dataset ```/resolutions/<resolution>/bins/<WEIGHTS_NAME>```.

9) ```"WEIGHTS_NAME": "weight"``` $\dots$ name of dataset in cooler containing normalization weights.

10) ```"SUB_M_SIZE_FIX": null``` $\dots$ fixed submatrix dimension.

11) ```"CHRSPLIT": 10``` $\dots$ number of submatrices into which the contact matrix is partitioned into.

12) ```"phi": 1``` $\dots$ number of bins to the next matrix.

13) ```"PHI_MAX": 1000``` $\dots$ number of submatrices; i.e. number of data points in entropy signal $S$. 
If set, $\varphi$ is increased until $\Phi \approx \Phi\_{\max}$.



# Running ENT3C

### Python:

* Command-Line Usage 
	* run ENT3C directly from terminal with: 

	```
	ENT3C <get_entropy|get_similarity|run_all> --config-file=/path/to/config_file/<config.json>
	```
	
	* ```<get_entropy>``` subcommand generate a dataframe with entropy values according to <config.json>. Output: ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_OUT.csv>```
	
	* ```<get_similarity>``` subcommand will generate a data frame with similarities according to <config.json> and ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_OUT.csv>```. Output: ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_similarity.csv```
	
	* ```<run_all>``` will generate both ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_OUT.csv>``` and ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_similarity.csv``` data frames. 

* or as python API 
		```
		import ENT3C
		ENT3C.run_get_entropy("config/config.json")
		ENT3C.run_get_similarity("config/config.json")
		ENT3C.run_all("config/config.json")
		```

### Julia:

* initial call for global package installation (see "initial julia set-up"): 
```
julia ENT3C.jl --config-file=config/config.test.json --install-deps=yes
```
* after that: 
```
julia ENT3C.jl --config-file=config/config.json
```

* alternativly load the predefined enviornments for julia 1.10.4 or 1.11.2  
```
julia ENT3C.jl --config-file=config/config.json --resolve-env=yes --julia-version=<v.v.v>
```

&#128161; note the matlab and julia implementations will always generate the entropy and similarity dataframes

### MATLAB:

```matlab -nodesktop -nosplash -nodisplay -r "ENT3C('config/config.json'); exit"```

Associated functions are contained in directories ```JULIA_functions/``` and ```MATLAB_functions/```.

&#128161; note the matlab and julia implementations will always generate the entropy and similarity dataframes

# Output files:
* ```<OUT_DIR>/<OUTPUT_PREFIX>_ENT3C_similarity.csv``` $\dots$ will contain all combinations of comparisons. The second two columns contain the short names specified in ```FILES``` and the third column ```Q``` the corresponding similarity score.  
	```
	Resolution	ChrNr	Sample1	Sample2	Q
	40000	2	HFFc6_BR2	A549_BR2	0.5584659814117208
	40000	2	HFFc6_BR2	G401_BR2	0.6594518933893059
	40000	2	HFFc6_BR2	HFFc6_BR1	0.8473530463515314
	.		.	.		.	.	.		.		.		.		.
	.		.	.		.	.	.		.		.		.		.
	.		.	.		.	.	.		.		.		.		.
	```

* ```<OUT_DIR>/<OUTPUT_PREFIX>_ENT3C_OUT.csv``` $\dots$ ENT3C output table. 
	```
	Name	ChrNr	Resolution	n	PHI	phi	binNrStart	binNrEND	START	END	S
	G401_BR1	2	40000	600	901	6	0	599	0	24000000	4.067424893091131
	G401_BR1	2	40000	600	901	6	6	605	240000	24240000	4.06198007393338
	G401_BR1	2	40000	600	901	6	12	611	480000	24480000	4.055473536905049
	.		.	.		.	.	.		.		.		.		.
	.		.	.		.	.	.		.		.		.		.
	.		.	.		.	.	.		.		.		.		.
	```
Each row corresponds to an evaluated submatrix with fields ```Name``` (the short name specified in ```FILES```), ```ChrNr```, ```Resolution```, the sub-matrix dimension ```sub_m_dim```, ```PHI=1+floor((N-SUB_M_SIZE)./phi)```, ```binNrStart``` and ```binNrEnd``` correspond to the start and end bin of the submatrix, ```START``` and ```END``` are the corresponding genomic coordinates and ```S``` is the computed von Neumann entropy.

```<OUT_DIR>/<OUTPUT_PREFIX>_ENT3C_signals.png``` $\dots$ simple visualization of entropy signals $S$:

Example entropy signals $S$ for ```ENT3C run_all --config=config/config.json``` (unbalanced 40kb contact matrices for even chromosomes across 5 cell lines):

<figure>
    <img src="OUTPUT/PYTHON/EvenChromosomes_NoWeights_40000_ENT3C_signals.png" style="max-width:80%;"
         alt="ENT3C python Output">
</figure>


# Data used in publication 
Both Julia and MATLAB implementations (```ENT3C.jl``` and ```ENT3C.m```) were tested on Hi-C and micro-C contact matrices binned at 40 kb in ```cool``` format. 

**micro-C** 
| Cell line | Biological Replicate (BR) | Accession (Experiemnt set) | Accession (```pairs```) |
| --- | --- | --- | --- |
| H1-hESC | 1 | 4DNES21D8SP8 | 4DNFING6ZFD, 4DNFIBMG8YA3, 4DNFIMT4PHZ1, 4DNFI8GM4EL9 |
| H1-hESC | 2 | 4DNES21D8SP8 | 4DNFIIYUGYBU, 4DNFI89L17XY, 4DNFIXP9MVBU, 4DNFI2YHYAJO, 4DNFIULY29IQ |
| HFFc6   | 1 | 4DNESphiT3UBH | 4DNFIN7IIIY6, 4DNFIJZDEIZ3, 4DNFIYBTHGNA, 4DNFIK8UIB5B |
| HFFc6   | 2 | 4DNESphiT3UBH | 4DNFIF5F4HRG, 4DNFIK82YRNM, 4DNFIATCW955, 4DNFIZU6ADT1, 4DNFIKWV6BY2  |
| HFFc6   | 3 | 4DNESphiT3UBH | 4DNFIFJL4JIH, 4DNFIONHB78N, 4DNFIG1ZOVIM, 4DNFIPKVL9YI, 4DNFIJM966UR, 4DNFIV8JNJB8 |

**Hi-C** 
| Cell line | Biological Replicate (BR) | Accession (Experiemnt set)  | Accession (```BAM```) |
| --- | --- | --- | --- |
| G401  | 1 | ENCSR079VIJ | ENCFF649MAY |
| G401  | 2 | ENCSR079VIJ | ENCFF758WUD |
| LNCaP | 1 | ENCSR346DCU | ENCFF977XHB |
| LNCaP | 2 | ENCSR346DCU | ENCFF204XII |
| A549  | 1 | ENCSR444WCZ | ENCFF867DCM |
| A549  | 2 | ENCSR444WCZ | ENCFF532XBC |

 1. for the Hi-C data, ```bam``` files were downloaded from the ENCODE data portal and converted into ```pairs``` files using the ```pairtools parse``` function<sup>3</sup>

	```pairtools parse --chroms-path hg38.fa.sizes -o <OUT.pairs.gz> --assembly hg38 --no-flip  --add-columns mapq  --drop-sam --drop-seq  --nproc-in 15 --nproc-out 15 <IN.bam>```

2. for the micro-C data, ```pairs``` of technical replicates (TRs) were merged with ```pairtools merge```. E.g. for H1-hESC, BR1 (4DNES21D8SP8):

	```pairtools merge -o <hESC.BR1.pairs.gz> --nproc 10 4DNFING6ZFDF.pairs.gz 4DNFIBMG8YA3.pairs.gz 4DNFIMT4PHZ1.pairs.gz 4DNFI8GM4EL9.pairs.gz```

3. 40 kb coolers were generated from the Hi-C/micro-C pairs files with ```cload pairs``` function<sup>4</sup>
 
	```cooler cload pairs -c1 2 -p1 3 -c2 4 -p2 5 --assembly hg38 <CHRSIZE_FILE:40000> <IN.pairs.gz> <OUT.cool>```

# References
1. Neumann, J. von., Thermodynamik quantenmechanischer Gesamtheiten. Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen. Mathematisch-Physikalische Klasse 1927. 1927. 273-291.
2. Felippe, H., et. al., Threshold-free estimation of entropy from a pearson matrix. EPL. 141(3):31003. 2023.
3. Open2C et. al., Pairtools: from sequencing data to chromosome contacts. bioRxiv. 2023. 
4. Abdennur,N., and Mirny, L.A., Cooler: scalable storage for Hi-C data and other genomically labeled arrays. Bioinformatics. 2020.
