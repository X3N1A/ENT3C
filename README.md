# ENT3C
Von Neuman entropy-based metric for 3C-Seq derived chromosomal contact matrices by comparing the "complexity" of patterns contained in local smaller submatrices along the diagonals of two contact matrices. 

## Requirements
Current implementation in MATLAB, Julia soon to follow.

Contact matrices in cool format (https://github.com/open2c/cooler)

! example.m works on 40 kb binned cool files !

# Data
example.m was tested on micro-C contact matrices of two biological replicates of the hESC (4DNFI9GMP2J8) cell-line and on Hi-C contact matrices of two biological replicates of the G401 (ENCSR079VIJ) (hg38)

 - download contact list-replicate (pairs) files for technical replicates corresponding to biological replicate
 - aggregate pairs using pairtools merge (https://github.com/open2c/pairtools)
 - generate 40kb cool files using cooler cload pairs
   
## Parameters
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
