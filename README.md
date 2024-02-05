# ENT3C
ENT3C is a method for qunatifying the similarity of 3C-Seq derived chromosomal contact matrices by comparing the "complexity" of patterns contained in smaller submatrices along their diagonals. It is based on the von Neumann entropy<sup>1</sup> and recent work for entropy quantification of Pearson correlation matrices<sup>2</sup>.

https://doi.org/10.1101/2024.01.30.577923 

## Summary of ENT3C approach
1. loads cooler files into MATLAB and looks for shared empty bins
2. ENT3C will extract smaller submatrices $\hat{a}$ of dimension $n\times n$ along the diagonal of an input contact matrix 
4. the logarithm of $\hat{a}$ is taken ($nan$s are set to zero)
5. $\hat{a}$ is transformed into a Pearson correlation matrix $\hat{P}$ ($nan$ values are set to zero)
6. $\hat{P}$ is transformed into $\hat{\rho}=\hat{P}/n$ to fulfill the conditions for computing the von Neumann entropy
7. the von Neumann entropy of $\hat{\rho}$ is computed as

   $S(\boldsymbol{\rho})=\sum_j \lambda_j \log \lambda_j$

   where $\lambda_j$ is the $j$ th eigenvalue of $\hat{\rho}$
8. this is repeated for subsequent submatrices along the diagonal of the input matrix and stored in the **"entorpy signal"** $S$
9. the Pearson correlation between $S$ of two matrices, is used as a similarity metric 

![explaination of ENT3C](Figures/ENT3C_explain.png)

**Figure:** Exemplary depiction of ENT3C derivation of the entropy signal $S$ of the contact matrix $\mathbf{A}$ of chromosome 14 binned at 40 kb of the HFFc6 cell line (biological replicate 1). ENT3C's parameters: submatrix dimension $n=300$, window shift $WS=10$, maximum number of data points in $S$, $WN_{MAX}=\infty$, were used, resulting in 146 submatrices along the diagonal of the contact matrix. For subsequent Pearson-transformed submatrices $\mathbf{P}$ along the diagonal of $\log{\mathbf{A}}$, ENT3C computes the von Neumann entropies $S_i(\mathbf{P_i})$; the resulting signal $S$ is shown in blue under the matrix. The first two ($\mathbf{P}_{1-2}$), the middle ($\mathbf{P}_{73}$), and the last two Pearson submatrices ($\mathbf{P}_{145-146}$) are shown.

# Requirements
Julia or MATLAB

# Data
Both Julia and MATLAB implementations (```ENT3C.jl``` and ```ENT3C.m```) were tested on Hi-C contact matrices in ```mcool```/```cool``` format of two biological replicates of the G401 (ENCSR079VIJ) and A549 (ENCSR444WCZ) cell-lines (hg38).

 1. download pairs files from ENCODE. G401 (ENCFF091BKE) A549 (ENCFF101MYU)
 2. generate cooler files with ```cload pairs``` function<sup>2</sup>
 3. generate multi-resolution mcool files with ```cload zoomify``` function<sup>2</sup> 
   
# Configuration Files
Both Julia and MATLAB implementations (```ENT3C.jl``` and ```ENT3C.m```) call configuration files in JSON format.

```config/config.julia.m```
```
{
  "WN_MAX": 1000,
  "CHRSPLIT": 7,
  "SUB_M_SIZE_FIX": null,
  "ChrNr": 14,
  "Resolution": 40000,
  "WS": 1,
  "NormM": 0,
  "DATA_PATH": "DATA_30e6",
  "FILES": [
    "ENCSR079VIJ.BioRep1.mcool",
    "G401_BR1",
    "ENCSR079VIJ.BioRep2.mcool",
    "G401_BR2",
    "ENCSR444WCZ.BioRep1.mcool",
    "A549_BR1",
    "ENCSR444WCZ.BioRep2.mcool",
    "A549_BR2"
  ],
  "OUT_DIR": "OUTPUT/JULIA",
  "OUT_PREFIX": "Chr14_40kb"
}
```
# Running main scripts 


# References
1. Neumann, J. von., Thermodynamik quantenmechanischer Gesamtheiten. Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse 1927, 1927, 273-291.
2. Felippe, H., et. al., Threshold-free estimation of entropy from a pearson matrix. EPL, 141(3):31003, 2023.
3. Abdennur,N., and Mirny, L.A., Cooler: scalable storage for Hi-C data and other genomically labeled arrays, Bioinformatics, 2020.
