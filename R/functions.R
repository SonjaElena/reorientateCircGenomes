library("devtools")
library(roxygen2)
library(dplyr)
library(tidyverse)
########################################
######### reorient fasta ###############
########################################
#' Reorientation of fasta files
#'
#' A function that batch reorientes a list of fasta files to start at the location given as vector of bp position or Protein IDs.
#' @param x list of DNAStringSets
#' @param replicon the replicon to be reoriented. Defaults to the largest replicon. Numbers increase with degreasing length (largest replicon = 1, second largest = 2, ...).
#' @param gff list of tables containing gff information. Required if location is given by Protein ID.
#' @param location vector of the starting positions in base pairs at which the Fasta file is to begin. Alternatively, a vector of protein IDs can be used if the corresponding gff objects are supplied.
#' @keywords genomics
#' @export
#' @examples
#' reorientFasta()

reorientFasta <- function(x, replicon = 1, gff = NA, location = position){

  if(is.na(gff) & is.character(location)) stop('gff file required if location is protein ID')

  if(!is.na(gff) == TRUE & is.character(location) == TRUE){
    location <- subset(gff, ProteinID == location)$start
    location <- do.call(rbind, location)
  }

  fastaReplicon <- x[width(x) == Rfast::nth(as.double(width(x)), replicon, descending = T)]

  info <- Seqinfo("Chromosome", seqlengths = width(fastaReplicon),
                  isCircular = T, genome = strsplit(names(fastaReplicon), split = " ")[[1]][1])


  fastaReoriented <- DNAStringSet(c(subseq(fastaReplicon, location,
                                           width(fastaReplicon))[[1]], subseq(fastaReplicon, 1, location-1)[[1]]))


  return(fastaReoriented)
}

#####################################
############ reorient gff ###########
#####################################
#' Reorientation of gff files
#'
#' A function that batch reorientes a list of gff files to start at the location given as vector of bp position or Protein IDs.
#' @param x list of tables containing gff information.
#' @param location starting positions in base pairs or of protein IDs.
#' @param replicon replicon identifiers based on which the Gff files should be subsetted.
#' @param prokka did the gff file originate as prokka output.
#' @keywords genomics
#' @export
#' @examples
#' reorientGff()

reorientGff <- function(x, location = position, replicon = NA, prokka = FALSE){

  if((sum(colnames(x) %in% c("V1"))<1) == TRUE) stop('Prokka must be set TRUE')

  if(is.na(replicon) ==  TRUE){
    chr_temp <- subset(x, V3 == "region")
    chr_temp <- chr_temp[which.max(chr_temp$V5), "V1"]
    gff_subset <- subset(x, V1 == chr_temp)
  } else {
    gff_subset <- x[x$V1 %in% replicon,]
  }

  chr_size <- rep(NA, length(gff_subset))
  chr_size <- subset(gff_subset, V3 == "region")$V5


  if(prokka == FALSE){
    gff_dataframe2 <- subset(gff_subset, V3 == "gene")
    gff_dataframe3 <- subset(gff_subset, V3 == "CDS")


    gff_split <- data.frame(ID = strsplit(gff_dataframe2$V9, split = ";") %>%
                              sapply(., function(x)
                                grep("^ID", x, value = T)) %>%
                              gsub("ID=", "", .),
                            LocusTag = strsplit(gff_dataframe2$V9, split = ";") %>%
                              sapply(., function(x)
                                grep("^locus_tag", x, value = T)) %>%
                              gsub("locus_tag=", "", .),
                            OldLocusTag = strsplit(gff_dataframe2$V9, split = ";") %>%
                              sapply(., function(x)
                                grep("old_locus_tag", x, value = T)) %>%
                              gsub("old_locus_tag=", "", .),
                            stringsAsFactors = FALSE)


    gff_split2 <- cbind(gff_dataframe2[,c(1,3,4,5,7)], gff_split)
    colnames(gff_split2)[1:5] <- c("Replicon", "Type", "start", "end", "strand")


    gff_splitB2<- data.frame(ID = strsplit(gff_dataframe3$V9, split = ";") %>%
                               sapply(., function(x) grep("Parent=", x, value = T)) %>%
                               gsub("Parent=", "", .),


                             ProteinID = strsplit(gff_dataframe3$V9, split = ";") %>%
                               sapply(., function(x) grep("protein_id=", x, value = T)) %>%
                               gsub("protein_id=", "", .),


                             Product = strsplit(gff_dataframe3$V9, split = ";") %>%
                               sapply(., function(x) grep("product", x, value = T)) %>%
                               gsub("product=", "", .),
                             stringsAsFactors = FALSE)

    gff <- inner_join(gff_splitB2, gff_split2, by = "ID")
  }else{gff <- gff_subset}


  if(is.character(location) == TRUE){
    location <- subset(gff, ProteinID == location)$start
  }

  gff$altend <- gff$end
  gff$altend[-nrow(gff)][gff$end[-nrow(gff)]>=gff$start[-1]] <- gff$start[-1][gff$end[-nrow(gff)]>=gff$start[-1]]-1

  gff$NewStart <- gff$start - (location-1)
  gff$NewEnd <- gff$altend - (location-1)

  testx <- subset(gff, gff$NewStart < 0)
  gff$NewStart[gff$NewStart < 0] <- chr_size+gff$NewStart[gff$NewStart < 0]

  testx <- subset(gff, gff$NewEnd < 0)
  gff$NewEnd[gff$NewEnd < 0] <- chr_size+gff$NewEnd[gff$NewEnd < 0]



  return(gff)
}

############################################
################## circular plot ###########
############################################
#' Circular genomic plots
#'
#' A function that plots the GC Skew, genes on plus and minus strand and location of selected genes. The output consists of a list containing the plot and a list of the regulators with their locations.
#' @param x list of DNAStringSets.
#' @param gff list of tables containing gff information. Files need to be output of function reorientGff() or colnames of start and end that should be used need to be labeled "NewStart" and "NewEnd".
#' @param proteinID vector of protein IDs whose gene position is to be shown in the plot.
#' @param radius radius adjustment of each circle.
#' @keywords genomics plot
#' @export
#' @examples
#' circGenomePlot()

circGenomePlot <- function(x, gff = NA, proteinID = proteinID, radius1 = 10, radius2 = 12, radius3 = 13, radius4 = 14, radius5 = 15){

  slw <- x %>% width(.) %>% seq(1,., 10000) %>% .[-length(.)] %>%
    IRanges(start = ., width = 10000)

  GC_Skew <- data.frame(start = start(slw),
                        end = end(slw),
                        skew = (letterFrequency(Views(x[[1]],
                                                      slw),
                                                letters = "G")-
                                  letterFrequency(Views(x[[1]],
                                                        slw),
                                                  letters = "C"))/
                          (letterFrequency(Views(x[[1]],
                                                 slw),
                                           letters = "G")+
                             letterFrequency(Views(x[[1]],
                                                   slw),
                                             letters = "C")))
  GC_Skew$type <- GC_Skew$G > 0

  genes_plus <- subset(gff, strand == "+")
  genes_minus <- subset(gff, strand == "-")

  regu <- gff[gff$ProteinID %in% proteinID,]

  info <- Seqinfo("Chromosome", seqlengths = max(GC_Skew$end) , isCircular = T, genome = "Rhodobacteraceae")

  gr_start <- GRanges(seqnames = "Chromosome",
                      IRanges(start = 0.01,
                              end = 1.01),
                      strand = "*",
                      seqinfo = info)

  gr_gcskew <- GRanges(seqnames = "Chromosome",
                       IRanges(start = GC_Skew$start,
                               end = GC_Skew$end+1),
                       strand = "*",
                       mcols = GC_Skew[,c(3,4)],
                       seqinfo = info)

  gr_plus <- GRanges(seqnames = "Chromosome",
                     IRanges(start = genes_plus$NewStart,
                             end = genes_plus$NewEnd),
                     strand = genes_plus$strand,
                     seqinfo = info)


  gr_minus <- GRanges(seqnames = "Chromosome",
                      IRanges(start = genes_minus$NewStart,
                              end = genes_minus$NewEnd),
                      strand = genes_minus$strand,
                      seqinfo = info)

  gr_reg <- GRanges(seqnames = "Chromosome",
                    IRanges(start = regu$NewStart,
                            end = regu$NewEnd),
                    strand = "*",
                    seqinfo = info_rueg)

  image <- ggbio() +
    circle(gr_plus,  geom = 'rect', space.skip = 0.0001, linetype = 0, fill = "steelblue",
           trackWidth = 1, radius = radius1) +

    circle(gr_minus,  geom = 'rect', space.skip = 0.0001, linetype = 0, fill = "steelblue",
           trackWidth = 1, radius = radius2) +

    circle(gr_start, fill =  "red", color =  "red",geom = "rect",
           radius = radius3, trackWidth = 1, linetype = 1, size = 0.5, space.skip = 0.0001) +

    circle(gr_reg, fill =  "black", color =  "black",geom = "rect",
           radius = radius4, trackWidth = 1, linetype = 1, size = 0.5, space.skip = 0.0001) +

    circle(gr_gcskew, geom = "bar", aes(fill = factor(mcols.type), y = mcols.G), size = 0.4,
           radius = radius5, trackWidth = 3, lty = "blank", space.skip = 0.0001) +

    scale_fill_manual(values = c("lightgrey", "#333333"), guide=FALSE)


  return(list(image, regu))
}


############################################
################## Gff tranform  ###########
############################################
#' Processing of Prokka Gff output
#'
#' A function that processes Gff files that originate from Prokka output.
#' @param x directory and file name of the Gff object to be processed. File could be read with function list.files.object.
#' @param reorient should function reorientGff() be used after this function = TRUE
#' @keywords Prokka gff
#' @export
#' @examples
#' gff_convert()


gff_convert <- function(x, reorient=FALSE){

  start_line <- readLines(x) %>% grep("##seq", .) %>% which.max() +1
  end_line <- readLines(x) %>% grep("##FASTA", .)

  gff_dataframe <- read.table(x,
                              nrows = end_line-start_line-1,
                              quote = "", sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE)

  gff_dataframe1 <- subset(gff_dataframe, V3 == "CDS")

  gff_split <- data.frame( Product = strsplit(gff_dataframe1$V9, split = ";") %>%
                             sapply(., function(x) grep("product", x, value = T)) %>%
                             gsub("product=", "", .),
                           LocusTag = strsplit(gff_dataframe1$V9, split = ";") %>%
                             sapply(., function(x)
                               grep("^locus_tag", x, value = T)) %>%
                             gsub("locus_tag=", "", .),
                           stringsAsFactors = FALSE)

  gff <- cbind(gff_dataframe1[,c(1,3,4,5,7)], gff_split)

  if(reorient == FALSE){
    colnames(gff)[1:5] <- c("Replicon", "Type", "start", "end", "strand")
  }

  rm(gff_dataframe, gff_dataframe1, gff_split, start_line, end_line)

  return(gff)

}
