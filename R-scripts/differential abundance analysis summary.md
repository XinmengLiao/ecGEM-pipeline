Totally 603 taxo, 46 samples, paired data (Day0 and Day 84)

### Number of significant species by different methods
| Method        | No. of P value < 0.05 | No. of P adj < 0.05 |
|---------------|-----------------------|---------------------|
| Wilcoxon (CLR) | 1                     | 0                   |
| DESeq2         | 362                   | 360                 |
| edgeR          | 119                   | 42                  |
| ANCOM-BC       | 85                    | 2                   |
| ALDEx2(we)     | 2                     | 1                   |
| ALDEx2(wi)     | 0                     | 0                   |
| MaAsLin2       | 6                     | 0                   |
**LefSe is not suitable for paired data. 


Totally 588 taxo (species levels) paired data (Day0 and Day 84)
**LefSe is not suitable for paired data. 
### Number of significant species by different methods (No filtration)
| Method        | No. of P value < 0.05 | No. of P adj < 0.05 | No. of sig increased metabolites | No. of sig decreased metabolites |
|---------------|-----------------------|---------------------|-----------------------|---------------------|
| Wilcoxon (CLR) | 1                     | 0                   |
| DESeq2         | 347                   | 347                 |176|103|
| edgeR          | 115                   | 42                  |198 (filtered out FDR = 1) | 133|
| ANCOM-BC       | 85                    | 2                   |
| ALDEx2(we)     | 2                     | 0                   |/|/|
| ALDEx2(wi)     | 0                     | 0                   |/|/|
| MaAsLin2       | 6                     | 0                   |/|/|


### Filtration 1: Exist both in Day0 and Day84
Totally 243 taxo (species levels) paired data (Day0 and Day 84)
| Method        | No. of P value < 0.05 | No. of P adj < 0.05 | No. of sig increased metabolites | No. of sig decreased metabolites |
|---------------|-----------------------|---------------------|-----------------------|---------------------|
| Wilcoxon (CLR) | 1                     | 0                   |
| DESeq2         | 35                   | 35                 |114 (filtered out FDR > 0.99) |85|
| edgeR          | 67                   | 38                  |161(filtered out FDR > 0.99)|135|
| ANCOM-BC       | 36                    | 1                   |

### Filtration 2: Exist in at least 5% of the samples
Totally 227 taxo (species levels) paired data (Day0 and Day 84)
| Method        | No. of P value < 0.05 | No. of P adj < 0.05 |
|---------------|-----------------------|---------------------|
| Wilcoxon (CLR) | 3                     | 0                   |
| DESeq2         | 72                   | 72                  |149(filtered out FDR > 0.99)|104|
| edgeR          | 96                   | 67                  |184 (filtered out FDR > 0.99)| 133|
| ANCOM-BC       | 46                    | 2                   |

### Results:
No filtration, filtered out FDR > 0.99: \
  for edgeR, only increased uptake network has sig L-/D-Serine, and L-/D-tartrate. Cystine has lower pvalue in increased than decreased. nicotinamide has lower pvalue in increased than decreased. But pvalues for cystine and nicotinamide are not sig. \
  for DESeq2, only increased network has sig L-Serine and D-Serine. But increased one has cystine which is not sig. L-/D-Tartrates are sig in increased only. 
