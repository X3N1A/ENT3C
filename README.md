# ENT3C
ENT3C is a method for qunatifying the similarity of 3C-Seq derived chromosomal contact matrices by comparing the "complexity" of patterns contained in smaller submatrices along their diagonals. It is based on the von Neumann entropy<sup>1</sup> and recent work for entropy quantification of Pearson correlation matrices<sup>2</sup>.

Current implementation of ENT3C is in MATLAB, Julia soon to follow. 

**summary of workflow**
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

## Requirements

Contact matrices in cool format (https://github.com/open2c/cooler)<sup>3</sup>

! example.m works on 40 kb binned cool files !

# Data
example.m was tested on micro-C contact matrices of two biological replicates of the hESC (4DNFI9GMP2J8) cell-line and on Hi-C contact matrices of two biological replicates of the G401 (ENCSR079VIJ) (hg38)

 - download contact list-replicate (pairs) files for technical replicates corresponding to biological replicate
 - aggregate pairs using pairtools merge (https://github.com/open2c/pairtools)<sup>4</sup>
 - generate 40kb cool files using cooler cload pairs
   
## Parameters ENT3C.m
```
input:
-CELL_TYPE ... name of cell-type to save in output table-Resolution ... resolution of cool/to extract from mcool file
-SUB_M_SIZE_FIX ... fixed submatrix size $n$
-CHRSPLIT ... determines window/sub-matrix size on which entropy values S(window) are calculated on
-WS ... shift size of submatrix alon diagonal
-MAX_WINDOWS maximum number of entropy values to compute (window shift is increased until desired window number is reached)
```
```
output:
-VN_ENT ... output table with Von Neuman entropy values of submatrices along diagonal and other information
```

## References
1. Neumann, J. von., Thermodynamik quantenmechanischer Gesamtheiten. Nachrichten von der Gesellschaft der Wissenschaften zu Göttingen, Mathematisch-Physikalische Klasse 1927, 1927, 273-291.
2. Felippe, H., et. al., Threshold-free estimation of entropy from a pearson matrix. EPL, 141(3):31003, 2023.
3. Abdennur,N., and Mirny, L.A., Cooler: scalable storage for Hi-C data and other genomically labeled arrays, Bioinformatics, 2020.
4. Open2C*, et. al., Pairtools: from sequencing data to chromosome contacts. bioRxiv, 2023.
