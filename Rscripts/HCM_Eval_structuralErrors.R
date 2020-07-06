#!/usr/bin/env Rscript
##############################################################
# Hybrid Correction Methods
# Evaluating structural errors
#------------------------------------------------------------#
# Date: April 2019
# Author: Lucus Pocus
##############################################################
# Overall statistics
# -> Nb of split reads (all spurious)
# -> Nb of accurate genomic alignments
# -> % of reads with at least one deleted exon [+ indication on del. exon length]
# -> % of reads with at least one inserted intron [+ indication on ins. exon length]
# Comparison to raw sequence
# -> Number of clarified exons
# -> Number of inserted errors that were not in raw sequence
##############################################################

str_folder <- "path/to/datasets"

aligners <- c("mmap2", "gmap", "graphmap")

ins_thr <- 20  #Minimal size for insertion in introns

##############################################################

suppressPackageStartupMessages( require(dplyr) )
suppressPackageStartupMessages( require(data.table) )
suppressPackageStartupMessages( require(stringr) )

##############################################################


##############################################################
# Overall statistics
#-------------------------------------#

# -> Nb of split reads (all spurious)
split.df <- data.frame()
# -> Nb of accurate genomic alignments
accurateAlign.df <- data.frame()
# -> % of Errors
Err <- data.frame()
# Comparison to raw sequence
# -> Number of clarified exons
# -> Number of inserted errors that were not in raw sequence

for(a in 1:length(aligners)){
  
  cat("Aligner:", aligners[a], "\n")
  ##################################
  
  file_pattern <- paste(".", aligners[a], ".struct_errors.txt", sep = "")
  
  files <- list.files(path = str_folder, full.names = TRUE, pattern = file_pattern)
  samples <- list.files(path = str_folder, full.names = FALSE, pattern = file_pattern)
  samples <- str_remove(samples, pattern = file_pattern)
  
for(f in 1:length(files)){
  
  print( files[f])
  dtmp <- data.frame( fread(file = files[f]) )
  
  #Spurious splits
  x <- dtmp %>% filter(split_read == TRUE) %>% distinct(read_name, .keep_all = TRUE) %>% nrow()
  split.df <- rbind.data.frame(split.df, data.frame(sample = samples[f], aligner = aligners[a], value = x))
  
  #Accurate gene-level alignments
  x1 <- dtmp %>% filter(split_read == FALSE & accurate_align == 1) %>% distinct(read_name, .keep_all = TRUE) %>% nrow()
  x2 <- dtmp %>% filter(split_read == FALSE) %>% distinct(read_name, .keep_all = TRUE) %>% nrow()
  accurateAlign.df <- rbind.data.frame(accurateAlign.df, 
                                       data.frame(sample = samples[f], aligner = aligners[a], nbAlign = x2, nbAccurate = x1))
  rm(x1,x2)
  
  nbEvalExons <- dtmp %>% filter(split_read == FALSE & accurate_align == 1 & interval_type == "exon") %>% nrow()
  nbEvalIntrons <- dtmp %>% filter(split_read == FALSE & accurate_align == 1 & interval_type == "introns") %>% nrow()
  
  x <- dtmp %>% filter(split_read == FALSE & accurate_align == 1) %>% 
                #filter(rank != 1 & rank != nbExons) %>% 
                group_by(read_name) %>% 
                summarize(nbDel = sum(inTx==1 & relative_overlap==0),
                          nbPerfectExon = sum(inTx==1 & relative_overlap==1),
                          nbImperfectExon = sum(inTx==1 & relative_overlap<1), 
                          nbIns = sum( as.numeric(interval_type == "intron" & relative_overlap*width>ins_thr)),
                          nbNoisy = sum( as.numeric(relative_overlap*width<=ins_thr & interval_type == "intron" & relative_overlap>0) )) %>% 
                data.frame()
  #Reads with del/ins/noisy things
  Err <- rbind.data.frame(Err, data.frame(sample = samples[f], 
                                          aligner = aligners[a], 
                                          nbEvalExons = nbEvalExons, 
                                          nbDel = sum(x$nbDel), nbReadWithDel = nrow(filter(x, nbDel>0)),
                                          nbPerfectExon = sum(x$nbPerfectExon),
                                          nbImperfectExon = sum(x$nbImperfectExon), 
                                          nbEvalIntrons = nbEvalIntrons,
                                          nbIns = sum(x$nbIns), nbReadIns = nrow(filter(x, nbIns>0)),
                                          nbReadwithErr = nrow(filter(x, nbDel>0 | nbIns>0)), 
                                          nbNoisy = sum(x$nbNoisy), nbReadNoise = nrow(filter(x, nbNoisy>0))))
  
 }
}

accurateAlign.df

Err2 <- Err %>% dplyr::select(sample, aligner, nbReadwithErr)
Err2 <- merge(accurateAlign.df, Err2, by = c("sample", "aligner"))
Err2 <-  mutate(Err2, proportion_error_free =  round((nbAccurate-nbReadwithErr)/nbAccurate,3))

View( Err2 )


##############################################################
# Comparison to raw data (taken as "control" data)
#-------------------------------------#

Err_ctrl <- data.frame()
for(a in 1:length(aligners)){
  
  cat("Aligner:", aligners[a], "\n")
  ##################################
  
  file_pattern <- paste(".", aligners[a], ".struct_errors.txt", sep = "")
  
  files <- list.files(path = str_folder, full.names = TRUE, pattern = file_pattern)
  samples <- list.files(path = str_folder, full.names = FALSE, pattern = file_pattern)
  samples <- str_remove(samples, pattern = file_pattern)
  
  control_file <- paste("/media/lucile/d72ec46a-1ede-4522-a59b-a8665ad566c2/TALC/Review/sim2_rand0/structure/sim2_s01.raw.", 
                        aligners[a], ".struct_errors.txt", sep = "")
  
  ifelse(length(grep(pattern = "raw", files))>0, treat_files <-  files[-grep(pattern = "raw", files)], treat_files <-  files)
  ifelse(length(grep(pattern = "raw", files))>0, treat_samples <- samples[-grep(pattern = "raw", files)], treat_samples <- samples)
  
  dtmp_ctrl <- data.frame( fread(file = control_file ) )

for(f in 1:length(treat_files)){
  
  print( treat_files[f] )
  dtmp <- data.frame( fread(file = treat_files[f]) )
  
  dtmp <- filter(dtmp,  split_read == FALSE & interval_type == "exon" & accurate_align==1)
  dtmp0 <- filter(dtmp_ctrl,  split_read == FALSE & interval_type == "exon" & accurate_align==1)
  
  dtmp <- filter(dtmp, read_name %in% dtmp0$read_name )
  dtmp0 <- filter(dtmp0, read_name %in% dtmp$read_name)
  
  dtmp <- arrange(dtmp, read_name, chr, start, end)
  dtmp0 <- arrange(dtmp0, read_name, chr, start, end)
  
  nrow(dtmp)
  nrow(dtmp0)
  
  stopifnot(length(setdiff(dtmp0$read_name, dtmp$read_name))== 0 )
  length(unique(dtmp0$read_name))
  length(unique(dtmp$read_name))
  
  dtmp$relative_overlap.ctrl <- dtmp0$relative_overlap
  
  x <- dtmp %>% filter(inTx==1 & relative_overlap.ctrl==1 & relative_overlap>0.9) %>% nrow()
  x0 <- dtmp %>% filter(inTx==1 & relative_overlap.ctrl==1) %>% nrow()
  y <- dtmp %>% filter(inTx ==1 & relative_overlap.ctrl==0 & relative_overlap>0.9) %>% nrow()
  y0 <- dtmp %>% filter(inTx ==1 & relative_overlap.ctrl==0) %>% nrow()
  Err_ctrl <- rbind.data.frame(Err_ctrl, 
                               data.frame(sample = treat_samples[f], 
                                          aligner = aligners[a],
                                          nbRawTrue = x0, 
                                          nbConserved=x, percentConserved=round(x/x0,3),
                                          nbRawFalse=y0,
                                          nbClarified = y, percentClarified=round(y/y0,3)))
 }
}

View( Err_ctrl )