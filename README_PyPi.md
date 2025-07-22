ENT3C is a method for qunatifying the similarity of micro-C/Hi-C derived chromosomal contact matrices. It is based on the von Neumann entropy<sup>1</sup> and recent work on entropy quantification of Pearson correlation matrices<sup>2</sup>.
For a contact matrix, ENT3C records the change in local pattern *complexity* of smaller Pearson-transformed submatrices along a matrix diagonal to generate a characteristic signal. Similarity is defined as the Pearson correlation between the respective entropy signals of two contact matrices.

https://github.com/X3N1A/ENT3C


## Requirements

* generate and activate python environment 
	
	```
	python3.12 -m venv .ent3c\_venv

	source .ent3c\_venv/bin/activate
	```

* install ENT3C and requirements via ```pyproject.toml```: 

	```
	pip install .
	```

  

## Running ENT3C

### Command-Line Usage 
* run ENT3C directly from terminal with: 

```
ENT3C <get_entropy|get_similarity|run_all> --config-file=/path/to/config_file/<config.json>
```
	
* ```<get_entropy>``` subcommand generate a dataframe with entropy values according to <config.json>. Output: ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_OUT.csv>```
	
* ```<get_similarity>``` subcommand will generate a data frame with similarities according to <config.json> and ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_OUT.csv>```. Output: ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_similarity.csv```
	
* ```<run_all>``` will generate both ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_OUT.csv>``` and ```OUTPUT/PYTHON/<OUT_PREFIX>_<_ENT3C_similarity.csv``` data frames. 

### or as python API 

```
import ENT3C
ENT3C.run_get_entropy("config/config.json")
ENT3C.run_get_similarity("config/config.json")
ENT3C.run_all("config/config.json")
```

## Parameters and configuration files of ENT3C

* The main ENT3C parameter affecting the final entropy signal $S$ is the dimension of the submatrices ```SUB_M_SIZE_FIX```. 

	* ```SUB_M_SIZE_FIX``` can be either be fixed by or alternatively, one can specify ```CHRSPLIT```; in this case ```SUB_M_SIZE_FIX``` will be computed internally to fit the number of desired times the contact matrix is to be paritioned into. 

	  ```PHI=1+floor((N-SUB_M_SIZE)./phi)```

	  where ```N``` is the size of the input contact matrix, ```phi``` is the window shift, ```PHI``` is the number of evaluated submatrices (consequently the number of data points in $S$).

* All implementations (```ENT3C.py```, ```ENT3C.jl``` and ```ENT3C.m```) use a configuration file in JSON format. 
	* example can be found in <config/config.json>

**ENT3C parameters defined in ```config/config.json```**
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

## Output files:
* ```<OUT_DIR>/<OUTPUT_PREFIX>_ENT3C_similarity.csv``` $\dots$ will contain all combinations of comparisons. The second two columns contain the short names specified in ```FILES``` and the third column ```Q``` the corresponding similarity score.  
```
Resolution	ChrNr	Sample1	Sample2	Q
cat OUTPUT/PYTHON/EvenChromosomes_NoWeights_ENT3C_similarity.csv  | head
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
cat OUTPUT/PYTHON/EvenChromosomes_NoWeights_ENT3C_similarity.csv  | head
Resolution	ChrNr	Sample1	Sample2	Q
Name	ChrNr	Resolution	n	PHI	phi	binNrStart	binNrEnd	START	END	S
G401_BR1	2	40000	600	901	6	0	599	0	24000000	4.067424893091131
G401_BR1	2	40000	600	901	6	6	605	240000	24240000	4.06198007393338
G401_BR1	2	40000	600	901	6	12	611	480000	24480000	4.055473536905049
G401_BR1	2	40000	600	901	6	18	617	720000	24720000	4.048004132456738
.		.	.		.	.	.		.		.		.		.
.		.	.		.	.	.		.		.		.		.
.		.	.		.	.	.		.		.		.		.
```
Each row corresponds to an evaluated submatrix with fields ```Name``` (the short name specified in ```FILES```), ```ChrNr```, ```Resolution```, the sub-matrix dimension ```sub_m_dim```, ```PHI=1+floor((N-SUB_M_SIZE)./phi)```, ```binNrStart``` and ```binNrEnd``` correspond to the start and end bin of the submatrix, ```START``` and ```END``` are the corresponding genomic coordinates and ```S``` is the computed von Neumann entropy.
