#' Processing of a gff file type downloaded from NCBI
#'
#' Process the gff file provided by NCBI into a data.frame.
#' @param x The unprocessed gff file downloaded from NCBI.
#' @return The processed gff file as data.frame with an additional column containing the alternative end, used for function(circGenomePlot)
#' @examples
#' gff_unprocessed <- path_to_gff
#' gff <- processNCBIgff(gff_unprocessed);
#' @export
processNCBIgff <- function(x){
  
  gff_dataframe <- read.table(x,
                              quote = "", sep = "\t",
                              header = FALSE,
                              stringsAsFactors = FALSE)
  
  gff_dataframe2 <- subset(gff_dataframe, V3 == "gene")
  gff_dataframe3 <- subset(gff_dataframe, V3 == "CDS")
  
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
  
  
  return(gff)
}


#' Processing of Prokka supplied gff file type
#'
#' Process the gff file created by Prokka into a data.frame.
#' @param x Path to the unprocessed gff file generated from Prokka.
#' @return The processed gff file as data.frame.
#' @examples
#' gff_unprocessed <- path_to_gff
#' gff <- processProkkagff(gff_unprocessed);
#' @export
processProkkagff <- function(x){
  
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
  
  colnames(gff)[1:5] <- c("Replicon", "Type", "start", "end", "strand")
  
  
  return(gff)
  
}





#' Reorientation of a gff file type
#'
#' Reorientation of the start position of a gff file type based on ProteinID or base pair position. New start and end locations are added in two additional columns called 'Ostart' and 'Oend'.
#' @param x gff file, processed with functions processProkkagff or processNCBIgff since the additional column called 'alternend' is used.
#' @param proteinID Supplies the start position at which file should be reoriented; Defaults to NA. If this ProteinID is not found on the biggest replicon, the replicon must be supplied as well.
#' @param bplocation Supplies the start position at which file should be reoriented in base pairs and needs to be identical with the start position of one of the proteins.
#' @param replicon Replicon to be reoriented. This option defaults to the largest replicon.
#' @param Rep_size Indicating the size of the Replicon to be used in base pairs. Alternatively, a genomic fasta sequence in fna format can be supplied.
#' @return The reoriented gff file with three additional columns called Ostart and Oend containing the adjusted start and end base pair locations as well as a column called 'Oaltend' containing the reoriented end position based on the alternative end column, in case the function 'processNCBIgff' has been used.
#' @examples
#' gff <- reorientgff(gff, "WP_012176686.1")
#' gff <- reorientgff(gff, bplocation = 1866, replicon = "CP000031.2")
#' gff <- reorientgff(gff, proteinID = "AAV97145.1", replicon = "CP000032.1")
#' gff <- reorientgff(gff, bplocation = 0, Rep_size = fna_path)
#' @export
reorientgff <- function(x, proteinID = NA, bplocation = bp_location, replicon = NA, Rep_size = fasta){
  
  
  # replicon
  if(is.na(replicon)){
    gff <- subset(x, Replicon == subset(data.frame(table(x$Replicon)), Freq == max(table(x$Replicon)))$Var1)
  }else{
    gff <- subset(x, Replicon == replicon)
  }
  
  # alternative end
  gff$alternend <- gff$end
  gff$alternend[-nrow(gff)][gff$end[-nrow(gff)]>=gff$start[-1]] <- gff$start[-1][gff$end[-nrow(gff)]>=gff$start[-1]]-1
  
  
  # location
  if(!is.na(proteinID)){
    bplocation <- subset(gff, ProteinID == proteinID)$start
  }
  
  
  # size of genome
  
  if(class(Rep_size) == "DNAStringSet"){
    names(Rep_size) <- names(Rep_size) %>% strsplit(., " ") %>% sapply(.,"[",1)
    chr_size2 <- Rep_size[names(Rep_size) %in% Rep_size[unique(gff$Replicon)]]
    chr_size <- width(chr_size2)
  }else{chr_size <- Rep_size}
  
  
  # reorientation
  gff$Ostart <- gff$start - bplocation
  gff$Oend <- gff$end - bplocation
  
  if(sum(gff$Ostart<0 | gff$Oend<0)>0){
    testx <- subset(gff, gff$Ostart < 0)
    gff$Ostart[gff$Ostart < 0] <- chr_size - (abs(testx$Ostart))
    
    testx <- subset(gff, gff$Oend < 0)
    gff$Oend[gff$Oend < 0] <- chr_size - (abs(testx$Oend))
  }
  
  
  # reorientation based on alternative end
  gff$Oaltend <- gff$alternend - bplocation
  testx <- subset(gff, gff$Oaltend < 0)
  gff$Oaltend[gff$Oaltend < 0] <- chr_size - (abs(testx$Oaltend))
  
  
  return(gff)
  
}


#' Reorientate a genome fasta file
#'
#' Adjust the start position of a fna file downloaded from NCBI.
#' @param x DNAStringSet of nucleotide sequences in fasta formate.
#' @param replicon Replicon to be reoriented. This option defaults to the largest replicon.
#' @param bplocation Location in base pairs that should be used as new start position.
#' @param proteinID ProteinID indicating the protein based on which the file should be reoriented. Must be accompanied by a processed gff file and either be located on the largest replicon or also be accompanied by an indication of the replicon to be used.
#' @param gff A processed gff file (e.g. using functions processProkkagff or processNCBIgff). Must be supplied when option ProteinID is selected.
#' @return
#' @examples
#' fasta <- reorientfna(x = dna_list, proteinID = "AAV93333.1", gff = gff3)
#' @export
reorientfna <- function(x, replicon = NA, bplocation = NA, proteinID = NA, gff = NA){
  
  # replicon
  if(is.na(replicon)){
    x <- x[width(x) %in% max(width(x))]
  }else{
    names(x) <- names(x) %>% strsplit(., " ") %>% sapply(.,"[",1)
    x <- x[names(x) %in% replicon]
  }
  
  # location
  if(!is.na(gff)){
    bp_location <- subset(gff, ProteinID == proteinID)$start
  }else{bp_location=bplocation}
  
  fasta_file <- DNAStringSet(c(unlist(subseq(x, bp_location, width(x))),
                               unlist(subseq(x, 1, bp_location-1))))
  
  
  return(fasta_file)
}


#' Generation of a circular genomic plot
#'
#' Generates a genomic plot indicating locations of regulators and showing the GC skew, based on gff file.
#' only one replicon.
#' @param x DNAStringSet of nucleotide sequences in fasta formate.
#' @param gff Processed gff file, has to be output of either functions processProkkagff() or processNCBIgff(), since a column with alternative end is generated that is used in this function. This file can have been reoriented afterwards with function reorientgff().
#' @param proteinID Vector of ProteinIDs to be indicated in the plot.
#' @param reorientgff If TRUE uses columns with names "Ostart" and "Oend" to obtain the base pair location. Defaults to FALSE and uses columns with names "start" and "end".
#' @return List containing the circular plot and data.frame of regulator location.
#' @examples
#' vect <- c("AAV93333.1", "AAV93335.1")
#' plot <- circGenomePlot(fasta, gff3, vect, reorigff = TRUE)
#' @export
circGenomePlot <- function(x, gff = gff, proteinID = proteinID, reorigff = FALSE){
  
  # GC skew
  slw <- x %>% width(.) %>% seq(1,., 5000) %>% .[-length(.)] %>%
    IRanges(start = ., width = 5000)
  
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
  
  
  
  
  # regulators
  genes_plus <- subset(gff, strand == "+")
  genes_minus <- subset(gff, strand == "-")
  
  regu <- gff[gff$ProteinID %in% proteinID,]
  
  # plot object
  info <- Seqinfo("Chromosome", seqlengths = max(GC_Skew$end) , isCircular = T, genome = "Genome")
  
  gr_start <- GRanges(seqnames = "Chromosome",
                      IRanges(start = 0.01,
                              end = 1.01),
                      strand = "*",
                      seqinfo = info)
  
  
  # GRange objects
  gr_gcskew <- GRanges(seqnames = "Chromosome",
                       IRanges(start = GC_Skew$start,
                               end = GC_Skew$end+1),
                       strand = "*",
                       mcols = GC_Skew[,c(3,4)],
                       seqinfo = info)
  
  if(reorigff == TRUE){
    gr_plus <- GRanges(seqnames = "Chromosome",
                       IRanges(start = genes_plus$Ostart,
                               end = genes_plus$Oaltend),
                       strand = genes_plus$strand,
                       seqinfo = info)
    gr_minus <- GRanges(seqnames = "Chromosome",
                        IRanges(start = genes_minus$Ostart,
                                end = genes_minus$Oaltend),
                        strand = genes_minus$strand,
                        seqinfo = info)
    gr_reg <- GRanges(seqnames = "Chromosome",
                      IRanges(start = regu$Ostart,
                              end = regu$Oaltend),
                      strand = "*",
                      seqinfo = info)
  }else{
    gr_plus <- GRanges(seqnames = "Chromosome",
                       IRanges(start = genes_plus$start,
                               end = genes_plus$alternend),
                       strand = genes_plus$strand,
                       seqinfo = info)
    gr_minus <- GRanges(seqnames = "Chromosome",
                        IRanges(start = genes_minus$start,
                                end = genes_minus$alternend),
                        strand = genes_minus$strand,
                        seqinfo = info)
    gr_reg <- GRanges(seqnames = "Chromosome",
                      IRanges(start = regu$start,
                              end = regu$alternend),
                      strand = "*",
                      seqinfo = info)
  }
  
  image <- ggbio() +
    circle(gr_plus,  geom = 'rect', space.skip = 0.0001, linetype = 0, fill = "steelblue",
           trackWidth = 1, radius = 9) +
    
    circle(gr_minus,  geom = 'rect', space.skip = 0.0001, linetype = 0, fill = "steelblue",
           trackWidth = 1, radius = 8) +
    
    circle(gr_start, fill =  "red", color =  "red",geom = "rect",
           radius = 10, trackWidth = 1, linetype = 1, size = 0.5, space.skip = 0.0001) +
    
    circle(gr_reg, fill =  "black", color =  "black",geom = "rect",
           radius = 10, trackWidth = 1, linetype = 1, size = 0.5, space.skip = 0.0001) +
    
    circle(gr_gcskew, geom = "bar", aes(fill = factor(mcols.type), y = mcols.G), size = 0.4,
           radius = 5, trackWidth = 3, lty = "blank", space.skip = 0.0001) +
    
    scale_fill_manual(values = c("lightgrey", "#333333"), guide=FALSE)
  
  
  return(list(image, regu))
}
