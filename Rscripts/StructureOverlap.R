#!/usr/bin/env Rscript
rm(list=ls())
args = commandArgs(TRUE)
##############################################################
# Hybrid Correction Methods
# Compare spliced alignments to reference exons
#------------------------------------------------------------#
# Date: April 2019
# Author: Lucus Pocus

##############################################################
# 
##############################################################
# Paths

#gtf_file <- "annotation.TxDb" # TxDb object build from gtf file
txdb_file <- "annotation.TxDb" # TxDb object build from gtf file

read_info_file <- "read_info.txt"

outFolder <- "/path/to/out"

if( length(args) < 2){
  stop("Two arguments must be supplied (input bam file and sample name) \n", call.=FALSE)
} else{
  bam_file <- args[1]
  sample <- args[2]
}

out_file <- paste(outFolder, "/", sample, ".struct_errors.txt", sep = "")

min_exonOverlap <- 5
min_intronOverlap <- 5

##############################################################
# R packages

suppressPackageStartupMessages( require(dplyr) )
suppressPackageStartupMessages( require(data.table) )
suppressPackageStartupMessages( require(stringr) )
suppressPackageStartupMessages( require(rtracklayer) )

suppressPackageStartupMessages( require(GenomicAlignments) )
suppressPackageStartupMessages( require(GenomicRanges) )
suppressPackageStartupMessages( library(GenomicFeatures) )

##############################################################

#txdb <- makeTxDbFromGFF(file = gtf_file)
#seqlevelsStyle(txdb) <- "UCSC"
#saveDb(txdb, file = txdb_file)

cat("[SV_eval]  Loading transcript database...")
txdb <- loadDb(txdb_file) 
seqlevelsStyle(txdb) <- "UCSC"

genes <- genes(x = txdb)
exons <- exonsBy(x = txdb, by = "tx", use.names = TRUE)
introns <- intronsByTranscript(x = txdb, use.names = TRUE)

cat("[SV_eval]  Loading sim. read description...")
read_info <- data.frame( fread(input = read_info_file) )
read_info$transcript_id <- tapply(read_info$trueRef_name, 
                                  INDEX = 1:nrow(read_info), 
                                  FUN = function(x) str_sub(string = x, end = str_length("ENST00000456328")))


cat("[SV_eval]  Loading alignments...\n ")
bam <- readGappedReads(file = bam_file, use.names = TRUE)
seqlevelsStyle(bam) <- "UCSC"

multi_map <- data.frame(read_name = names(bam)) %>% group_by(read_name) %>% summarize(nb_map = n()) %>% data.frame()

cat("[SV_eval]  Stats for split alignments: \n")
print( table(multi_map$nb_map) )

# read_name transcript_id gene_id nbIns nbDel
result.df <- data.frame()

cat("[SV_eval] Checking input alignments...\n")
for(read in 1:length(bam)){
  read_name <- names( bam[ read ] )
  transcript_id <- read_info$transcript_id[ which(read_info$read_name == read_name)]
  gene_id <- read_info$gene_id[ which(read_info$read_name == read_name)]
  gene_id <- str_sub(string = gene_id, end = str_length("ENSG00000000003"))
  
  nb_align <-  multi_map$nb_map[ which(multi_map$read_name == read_name) ]
  accurate_align <- FALSE
  
  df.read <- data.frame()
  
  valid <- (gene_id %in% names(genes))
  if( !valid ) cat("Invalid gene found:", gene_id, "\n")
  valid <- ( valid & (transcript_id %in% names(exons)) )
  if( !valid ) cat("Invalid transcript found:", transcript_id, "\n")
  
  if( valid ){
    
    chr <- unique(as.vector(seqnames(exons[ transcript_id ])[[1]]))
    pex <- sort(unlist(ranges( exons[ transcript_id ], use.mcols = TRUE) ))
    pex.df <- data.frame( mcols(pex) )
    
    pint <- unlist(ranges( introns[ transcript_id ], use.mcols = TRUE))
    
    #First check that read overlaps with reference gene:
    accurate_align <- NA
    hit <- suppressWarnings( findOverlaps( bam[ read ], genes[ gene_id ]) )
    accurate_align <- (length(hit)>0)
    
    #spliced alignment
    LR_ranges <- extractAlignmentRangesOnReference(cigar = cigar( bam[ read ] ), 
                                                   pos = start( bam[ read ] ), 
                                                   drop.D.ranges = FALSE, f = NULL)
    LR_ranges <- LR_ranges[[1]]; LR_ranges
    
    hits.pex <- findOverlaps(query = LR_ranges, subject = pex, minoverlap = min_exonOverlap)
    pex.df$hit <- 0
    pex.df$hit[ unique(subjectHits(hits.pex)) ] <- 1
    
    pex.df$relative_overlap <- 0
    relative_overlap <- width( pintersect(LR_ranges[ queryHits(hits.pex)], pex[ subjectHits(hits.pex) ]) ) / width(pex[ subjectHits(hits.pex) ] )
    relative_overlap <- round(relative_overlap, 2)
    pex.df$relative_overlap[ subjectHits(hits.pex) ] <- relative_overlap
    
    df.read <- data.frame(read_name = read_name, gene_id = gene_id, transcript_id = transcript_id, 
                          accurate_align = as.numeric(accurate_align),
                          chr = chr, start = start(pex), end = end(pex), width = width(pex),
                          interval_type = "exon", rank = pex.df$exon_rank, inTx = 1, inReadAlign = pex.df$hit, relative_overlap = pex.df$relative_overlap)
    
    if( length(pint)>0 ){
      pint.df <- data.frame( mcols(pint) )
      hits.pint<- findOverlaps(query = LR_ranges, subject = pint, minoverlap = min_exonOverlap)
      pint.df$hit <- 0
      pint.df$relative_overlap <- 0
      if( length(hits.pint)>0 ){
        pint.df$hit[unique(subjectHits(hits.pint))] <- 1
        
        relative_overlap <- width( pintersect(LR_ranges[queryHits(hits.pint)], pint[ subjectHits(hits.pint) ]) ) / width(pint[ subjectHits(hits.pint) ] )
        relative_overlap <- round(relative_overlap, 5)
        pint.df$relative_overlap[ subjectHits(hits.pint) ] <- relative_overlap
      }
      df.read <- rbind.data.frame(df.read,  
                                  data.frame(read_name = read_name, gene_id = gene_id, transcript_id = transcript_id, 
                                             accurate_align = as.numeric(accurate_align),
                                             chr = chr, start = start(pint), end = end(pint), width = width(pint), interval_type = "intron", 
                                             rank = 1:nrow(pint.df), inTx = 0, inReadAlign = pint.df$hit, relative_overlap = pint.df$relative_overlap))
    }
    df.read$nbExons <- nrow(pex.df)
    df.read$split_read <- (nb_align>1)
    #df.read <- arrange(df.read, start, end)
    result.df <- rbind.data.frame(result.df, df.read)
  } #END-IF not valid reference names
} #END FOR on reads

fwrite(x = result.df, file = out_file)


