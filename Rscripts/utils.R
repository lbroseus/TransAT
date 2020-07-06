##############################################################
# Functions used in Evaluation scripts 
#------------------------------------------------------------#
# Copyright (C) 2019-2020   <lucile.broseus@igh.cnrs.fr>
##############################################################

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(GenomicAlignments)
  library(dplyr)
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(compositions)
})


#_______________________________________________________#
#' computeBaseAccuracy() 
#' 
#' @description
#' This functions computes basic alignment statitics (#/% of alignments)
#' and estimates base accuracy from primary alignments as:
#' \[ \dfrac{nbrM+nbrI+nbrD−NM}{nbrM+nbrI+nbrD} \] 
#' 
#' @param genomebamdir Folder containing alignments in bam files
#' @param file_pattern Specific pattern given to the bam files
#'                          
#' @return A data.frame             
#--------------------------------------------------------#
# Partly copied from:
# https://github.com/csoneson/NativeRNAseqComplexTranscriptome/tree/master/Rscripts

computeBaseAccuracy <- function(genomebamdir, file_pattern){
  
  bamfiles <- list.files(genomebamdir, pattern = file_pattern, recursive = F, full.names = TRUE)
  
  nbPrimaryAlignments <- 0
  nbReads <- 0

  #cat("Detected bam files to be analysed: \n"); print( bamfiles )
  
  samples <- basename(bamfiles) %>% stringr::str_remove(pattern = file_pattern)
  
  summaries <- matrix(nrow = length(bamfiles), ncol = 11)
  
  #Compute base accuracy as described in Soneson2019 for all bam files
  for(f in 1:length(bamfiles)){
    
    nbReads <- countBam(bamfiles[f], param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = TRUE)))$records
    
    #cat(bamfiles[f], "\n")
    bam <- readGAlignments(bamfiles[f], use.names = TRUE, 
                           param = ScanBamParam(tag = c("NM"),
                                                what = c("qname","flag", "rname", "pos", "mapq")))
    
    #Keep only aligned reads and their primary alignment
    bam <- subset(bam, flag %in% c(0, 16))
    
    nbPrimaryAlignments <- length( unique(names(bam)) )
    nbReads <- nbReads + nbPrimaryAlignments
    
    #Redundant, but at least, it works...
    
    opts <- "S"
    widths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = opts)
    mcols(bam)$nbrS <- sapply(widths, sum)
    
    opts <- "H"
    widths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = opts)
    mcols(bam)$nbrH <- sapply(widths, sum)
    
    opts <- "M"
    widths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = opts)
    mcols(bam)$nbrM <- sapply(widths, sum)
    
    opts <- "I"
    widths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = opts)
    mcols(bam)$nbrI <- sapply(widths, sum)
    
    opts <- "D"
    widths <- GenomicAlignments::explodeCigarOpLengths(cigar(bam), ops = opts)
    mcols(bam)$nbrD <- sapply(widths, sum)
    
    mcols(bam)$readLength <- rowSums( as.matrix(mcols(bam)[, c("nbrS", "nbrH", "nbrM", "nbrI")]) )
    
    mcols(bam)$base_accuracy <- (mcols(bam)$nbrM+mcols(bam)$nbrI+mcols(bam)$nbrD-mcols(bam)$NM)/(mcols(bam)$nbrM+mcols(bam)$nbrI+mcols(bam)$nbrD)
    
    mcols(bam)$Mis <- ifelse(mcols(bam)$NM==0, 0, (mcols(bam)$NM-mcols(bam)$nbrI-mcols(bam)$nbrD)/(mcols(bam)$NM)*(1-mcols(bam)$base_accuracy))
    mcols(bam)$Ins <- ifelse(mcols(bam)$NM==0, 0,(mcols(bam)$nbrI)/(mcols(bam)$NM)*(1-mcols(bam)$base_accuracy))
    mcols(bam)$Del <- ifelse(mcols(bam)$NM==0, 0,(mcols(bam)$nbrD)/(mcols(bam)$NM)*(1-mcols(bam)$base_accuracy))
    
    
    summaries[f, ] <- c(nbPrimaryAlignments, 
                        percentAlignments = round(nbPrimaryAlignments/nbReads*100,4),
                        summary( round(mcols(bam)$base_accuracy*100,4) ),
                        round(mean(mcols(bam)$Mis*100), 4), 
                        round(mean(mcols(bam)$Ins*100), 4),
                        round(mean(mcols(bam)$Del*100), 4))
  
  }
  
  colnames(summaries) <- c("nbPrimaryAlignments", 
                           "PercentAlignments",
                           "Min", "1stQ", "Median", "Mean", "3rdQ", "Max",
                           "Mis", "Ins", "Del")
  
  summaries <- data.frame(sample = samples, summaries)
  
  return( summaries )
  
}
