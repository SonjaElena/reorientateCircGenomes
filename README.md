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

``` C
x	
```
list of DNAStringSets

replicon = the replicon to be reoriented. Defaults to the largest replicon. Numbers increase with degreasing length (largest replicon = 1, second largest = 2, ...).

gff	= list of tables containing gff information. Required if location is given by Protein ID.

location	= vector of the starting positions in base pairs at which the Fasta file is to begin. Alternatively, a vector of protein IDs can be used if the corresponding gff objects are supplied.


``` C
reorientGff(x, location = position, replicon = NA, prokka = FALSE)
```

``` C
circGenomePlot <- function(x, gff = NA, proteinID = proteinID, radius1 = 10, radius2 = 12, radius3 = 13, radius4 = 14, radius5 = 15)
```

``` C
gff_convert <- function(x, reorient=FALSE)
```
