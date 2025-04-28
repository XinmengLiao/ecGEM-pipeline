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
### Results (p-corrected in uptake network)
| Method           | DA                         | D-serine | L-serine | Cysteine | D-tartrate | L-tartrate |
|------------------|-----------------------------|----------|----------|---------|------------|------------|
| rm.org.up        | edgeR                       | 0.25     | 0.036 *  | 0.013 * | 0.83       | 0.50       |
| rm.org.down      | edgeR                       | 0.80     | 0.079    | 0.74    | 0.99       | 0.98       |
| rm.org.up        | edgeR (filtered pvalue<1)    | 0.009    | 0.13     | 0.76    | 0.59       | 0.12       |
| rm.org.down      | edgeR (filtered pvalue<1)    | 0.46     | 0.92     | 0.37    | 0.99       | 0.96       |
| rm.org.up        | DESeq2                      | 0.047 *  | 0.007 *  | 0.27    | 0.24       | 0.026      |
| rm.org.down      | DESeq2                      | 0.053    | 0.32     | 0.96    | 0.93       | 0.90       |
| rm.org.up        | DESeq2 (filtered pvalue<0.99)| 0.051    | 0.005 *  | 0.32    | 0.21       | 0.024 *    |
| rm.org.down      | DESeq2 (filtered pvalue<0.99)| 0.06     | 0.31     | 0.90    | 0.35       | 0.91       |


