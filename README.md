# reorient_circular_Genomes
![image-01](https://user-images.githubusercontent.com/51213363/89195191-f9b3db00-d582-11ea-9638-7cbf209d9162.png)
reorient_circular_Genomes is a R package to modify and visualize circular genomes. It can reorient fasta and gff files that originate from NCBI or Prokka based on base pair location or ProteinID. It can also visualize the circular genome, including strand information, the GC Skew and locations of selected genes.

This package has first been introduced here:
[**Koppenhöfer S, Tomasch J, Lang A.S.** Title. *Journal* year.](link here)

## Dependencies
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [BSgenome](http://bioconductor.org/packages/release/bioc/html/BSgenome.html)
- [tidyverse](https://www.tidyverse.org/)
- [ggbio](http://www.bioconductor.org/packages/release/bioc/html/ggbio.html/)


## Installation

## Usage
``` C
reorientFasta(x, replicon = 1, gff = NA, location = position)
```
