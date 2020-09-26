# reorient_circular_Genomes
![image-01](https://user-images.githubusercontent.com/51213363/89195191-f9b3db00-d582-11ea-9638-7cbf209d9162.png)
reorient_circular_Genomes is a R package to modify and visualize circular genomes. It can reorient fasta and gff files that originate from NCBI or Prokka based on base pair location or ProteinID. It can also visualize the circular genome, including strand information, the GC Skew and locations of selected genes.

This package has first been introduced here:
[**Koppenhöfer S., Tomasch J., Lang A.S.** Title. *Journal* year.](link here)

## Dependencies
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
- [BSgenome](http://bioconductor.org/packages/release/bioc/html/BSgenome.html)
- [tidyverse](https://www.tidyverse.org/)
- [ggbio](http://www.bioconductor.org/packages/release/bioc/html/ggbio.html/)


## Installation

## Usage

A function that batch reorientes a list of fasta files to start at the location given as vector of bp position or Protein IDs.
``` C
reorientFasta(x, replicon = 1, gff = NA, location = position)
```
* `x` a list of DNAStringSets

* `replicon` the replicon to be reoriented. Defaults to the largest replicon. Numbers increase with degreasing length (largest replicon = 1, second largest = 2, ...).

* `gff`	a list of tables containing gff information. Required if location is given by Protein ID.

* `location` vector of the starting positions in base pairs at which the Fasta file is to begin. Alternatively, a vector of protein IDs can be used if the corresponding gff objects are supplied.


A function that batch reorientes a list of gff files to start at the location given as vector of bp position or Protein IDs.
``` C
reorientGff(x, location = position, replicon = NA, prokka = FALSE)
```
* `x` a list of tables containing gff information.

* `location` starting positions in base pairs or of protein IDs.

* `replicon` replicon identifiers based on which the Gff files should be subsetted.

* `prokka` did the gff file originate as prokka output.


A function that plots the GC Skew, genes on plus and minus strand and location of selected genes. The output consists of a list containing the plot and a list of the regulators with their locations.
``` C
circGenomePlot(x, gff = NA, proteinID = proteinID, radius1 = 10, radius2 = 12, radius3 = 13, radius4 = 14, radius5 = 15)
```
* `x` a list of DNAStringSets.
* `gff` list of tables containing gff information. Files need to be output of function reorientGff() or colnames of start and end that should be used need to be labeled "NewStart" and "NewEnd".
* `proteinID` vector of protein IDs whose gene position is to be shown in the plot.
* `radius` allows the adjustment of each circle.


A function that processes Gff files that originate from Prokka output.
``` C
gff_convert(x, reorient=FALSE)
```
* `x` directory and file name of the Gff object to be processed. File could be read with function list.files.object.
* `reorient` should function reorientGff() be used after this function = TRUE

