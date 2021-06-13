library(ShortRead)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(IRanges)

tmpFile <- function(){ paste0(paste0(stringi::stri_rand_strings(30, 1, '[A-Za-z0-9]'), collapse = ''), '.tmp') }


# Read in sample metadata.
samples <- readr::read_delim('sampleConfig.inward.allSamples.tsv', '\t')


# Configuration.
R2.analysis.width <- 200
R2.recognitionSeqWidth <- 32
minQualityScore <- '5'
blastn <- '/home/everett/blast+/bin/blastn'
cutadapt <- '/usr/bin/cutadapt3'
p5.expected.200 <- 'CATGCTCTAGGAAGATCCAGGTTAATTTTTAAAAAGCAGTCAAAAGTCCAAGTGGCCCTTGGCAGCATTTACTCTCTCTGTTTGCTCTGGTTAATAATCTCAGGAGCACAAACATTCCAGATCCAGGTTAATTTTTAAAAAGCAGTCAAAAGTCCAAGTGGCCCTTGGCAGCATTTACTCTCTCTGTTTGCTCTGGTTAA'
p3.expected.200 <- 'CATGCTCTAGGAAGATCCTTATCGATTTTACCACATTTGTAGAGGTTTTACTTGCTTTAAAAAACCTCCCACATCTCCCCCTGAACCTGAAACATAAAATGAATGCAATTGTTGTTGTTAACTTGTTTATTGCAGCTTATAATGGTTACAAATAAAGCAATAGCATCACAAATTTCACAAATAAAGCATTTTTTTCACTG'
invisible(file.remove(list.files('./tmp', full.names = TRUE)))


# Read in seq data.
I1 <- readFastq('AAV_in_I1.fastq.gz')
R1 <- readFastq('AAV_in_R1.fastq.gz')
R2 <- readFastq('AAV_in_R2.fastq.gz')


# Remove lane info from read ids to make them comparable. 
I1@id <- BStringSet(sub('\\s.+$', '', as.character(I1@id)))
R1@id <- BStringSet(sub('\\s.+$', '', as.character(R1@id)))
R2@id <- BStringSet(sub('\\s.+$', '', as.character(R2@id)))


# Trim common linker from over-read fragments.
# This is done before quality trimming in order to preserve possibly low quality adapter sequences.
suppressWarnings(file.remove(c('R2.ff', 'R2.ff.out', 'R2.log')))
writeFastq(R2, file = 'R2.ff', compress = FALSE)
system(paste0(cutadapt, ' -f fastq  -e 0.15 -a AGTCCCTTAAGCGGAG --overlap 2 --info-file=R2.log R2.ff > R2.ff.out'))
R2 <- readFastq('R2.ff.out')
suppressWarnings(file.remove(c('R2.ff', 'R2.ff.out')))


# Quality trim R2 reads.
R2 <- trimTailw(R2, 2, minQualityScore, 5)
R2 <- R2[width(R2) >= R2.analysis.width]


# Sync reads since a number of R2 reads were likely removed.
I1 <- I1[I1@id %in% R2@id]
R1 <- R1[R1@id %in% R2@id]
if(! all(R2@id == I1@id)) stop('Error - R2 and I1 reads are not synced.')


# Cycle through reach row in the sample table, demultiplex reads, identify sequencing direction,
# and return subseq of reads.
r <- bind_rows(lapply(split(samples, 1:nrow(samples)), function(x){
  # Demultiplex reads for this sample table row using both unique linker sequences and I1 barcodes. 
  a <- as.character(subseq(R1@sread, 1, 20)) == substr(x$adriftRead.linker.seq, 1, 20)
  b <- as.character(I1@sread) == as.character(reverseComplement(DNAString(x$index1.seq)))
  reads <- R2[a & b]
  
  if(length(reads) == 0) return(data.frame())
  
  # Define the first 50 NTs of expected 5' and 3' sequences as recognition sequences so we know
  # which reaction (5' or 3') we are seeing.
  p5.leader <- substr(p5.expected.200, 1, R2.recognitionSeqWidth)
  p3.leader <- substr(p3.expected.200, 1, R2.recognitionSeqWidth)
  
  reads.5p <- reads[vcountPattern(p5.leader, subseq(reads@sread, 1, R2.recognitionSeqWidth), max.mismatch = round(R2.recognitionSeqWidth/10)) > 0]
  reads.3p <- reads[vcountPattern(p3.leader, subseq(reads@sread, 1, R2.recognitionSeqWidth), max.mismatch = round(R2.recognitionSeqWidth/10)) > 0]
  
  tibble(subject = x$subject, sample = x$sample, replicate = x$replicate,
         reads_5p = list(reads.5p), reads_3p = list(reads.3p))
}))


blastReads <- function(reads, db){
  f <- tmpFile()
  writeFasta(reads, paste0('tmp/', f))
  
  system(paste0(blastn, ' -num_threads 30 -word_size 7 -evalue 50 -outfmt 6 -query tmp/', f, ' -db ', db, ' -out tmp/', f, '.blast'))
  b <- data.frame()
  if(file.info(paste0('tmp/', f, '.blast'))$size > 0){
    b <- read.table(paste0('tmp/', f, '.blast'), sep = '\t', header = FALSE)
    names(b) <- c('qname', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')
    b$alignmentLength <- b$qend - b$qstart + 1
    b$pcoverage <- (b$alignmentLength / R2.analysis.width) * 100
  }     
  invisible(file.remove(list.files('tmp', pattern = f, full.names = TRUE)))
  return(b)
}

blast2rearangements <- function(x){
  if(nrow(x) == 0) return(data.frame())
  x <- subset(x, alignmentLength >= 15 & evalue <= 1e-10)
  if(nrow(x) == 0) return(data.frame())
  
  library(parallel)
  n <- 22
  cluster <- makeCluster(n)
  
  # Here we create a spliting variable across the blast data frame
  # being careful not to split read ids into different chunks.
  a <- floor(n_distinct(x$qname) / n)
  b <- 1
  
  o <- bind_rows(lapply(split(x, x$qname), function(x2){
         x2$n <- b
         b <<- b + 1
         if(b > a) b <<- 1
         x2
       }))
  
  r <- bind_rows(parLapply(cluster, split(o, o$n), function(b){
         library(dplyr)
         library(IRanges)
    
         bind_rows(lapply(split(b, b$qname), function(b2){
           b2 <- arrange(b2, qstart)
           b2$strand <- ifelse(b2$send < b2$sstart, '-', '+')
      
           b2$sstart2 <- ifelse(b2$sstart > b2$send, b2$send, b2$sstart)
           b2$send2   <- ifelse(b2$sstart > b2$send, b2$sstart, b2$send)
      
           # Shrink the ranges to prevent closely spaced ranges from being assembeled.
          b2$qstart <- b2$qstart + 3
          b2$qend   <- b2$qend - 3
      
          ir <- IRanges(start = b2$qstart, end = b2$qend)
          if(length(ir) == 0) return(data.frame())
      
          names(ir) <- paste0(b2$qstart, '..', b2$qend, '[', b2$sstart2, b2$strand, b2$send2, ']')
      
          o <- ir[1]
          invisible(lapply(split(ir, 1:length(ir)), function(a){
            if(all(! countOverlaps(o, a) > 0)){
            o <<- c(o, a)
          }
         }))
      
         if(length(o) == 0) return(data.frame())
      
         # Undo the range shrinkage.
         n1 <- as.integer(stringr::str_extract(unlist(lapply(strsplit(names(o), '\\['), '[', 1)), '^\\d+')) - 3
         n2 <- as.integer(stringr::str_extract(unlist(lapply(strsplit(names(o), '\\['), '[', 1)), '\\d+$')) + 3
         names(o) <- paste0(n1, '..', n2, stringr::str_extract(names(o), '\\[.+\\]'))
      
         data.frame(readID = b2$qname[1], rearrangement = paste0(names(o), collapse = ';'))
       }))
    }))
  
    stopCluster(cluster)
    r
}


readTable <- function(x){
   tibble(id = as.character(x@id), readSeq = as.character(x@sread)) %>% 
   dplyr::group_by(readSeq) %>%
   dplyr::summarise(readID = id[1], n = n_distinct(id)) %>%
   dplyr::ungroup() %>%
   dplyr::select(-readSeq)
}

# Collapse replicate reads.
r2 <- lapply(split(r, r$subject), function(x){
  message(x$sample[1])   
  
  p5 <- Reduce('append', x$reads_5p)
  p3 <- Reduce('append', x$reads_3p)
  
  if(length(p5) == 0 | length(p3) == 0) return(data.frame())
  
  if(grepl('TBG-HC', x$subject[1])){
    db <- 'sequences/heavyChainVector'
  } else if(grepl('TBG-LC', x$subject[1])){
    db <- 'sequences/lightChainVector'
  } else {
    db <- 'sequences/vectorSequences'
  }
  
  
  p5.tbl <- readTable(p5)
  p5 <- p5[as.character(p5@id) %in% p5.tbl$readID]
  r.p5 <-blast2rearangements(blastReads(p5, db))
  message('p5 rearangements done')   
  
  p3.tbl <- readTable(p3)
  p3 <- p3[as.character(p3@id) %in% p3.tbl$readID]
  r.p3 <-blast2rearangements(blastReads(p3, db))
  message('p3 rearangements done')  
  
  d <- bind_rows(p5.tbl, p3.tbl)
  r <- bind_rows(r.p5, r.p3)
  
  r$subject <- x$subject[1]
  r$sample <- x$sample[1]
  
  list(d, r)
})


tbls <- bind_rows(lapply(r2, function(x){ 
         if(length(x) == 2){
           return(x[[1]])
         } else {
           return(tibble())
         }}))

rcbs <- bind_rows(lapply(r2, function(x){ 
          if(length(x) == 2){
            return(x[[2]])
          } else {
            return(tibble())
         }})) %>% left_join(tbls, by = 'readID') %>%
         group_by(subject, rearrangement) %>%
         summarise(reads = sum(n)) %>%
         ungroup()

rcbs <- bind_rows(lapply(split(rcbs, 1:nrow(rcbs)), function(x){
          o <- unlist(strsplit(x$rearrangement, ';'))
          x$length <- as.integer(str_extract(str_extract(o[length(o)], '\\.\\.\\d+'), '\\d+'))
          x$n <- length(o)
          x
        })) %>%
        dplyr::filter(length >= 200)

readCounts <- group_by(rcbs, subject) %>%
              summarise(totalReads = sum(reads)) %>%
              ungroup()

rcbs2 <- group_by(dplyr::filter(rcbs, n > 1 & reads > 1), subject, rearrangement) %>%
         summarise(totalReads = sum(reads)) %>%
         ungroup()

rcbs2 <- bind_rows(lapply(split(rcbs2, 1:nrow(rcbs2)), function(x){
          o <- unlist(strsplit(x$rearrangement, ';'))
          a <- as.integer(unlist(str_extract_all(x$rearrangement, '\\d+')))
          x$a <- a[2]
          x$b <- a[5]
          a <- unlist(str_extract_all(x$rearrangement, '[+-]'))
          x$as <- a[1]
          x$bs <- a[2]
          x
        })) %>% arrange(subject, a, as, b, bs) %>% select(-a, -as, -b, -bs)

write.table(rcbs2, file = 'rcbs', sep = '\t', row.names = FALSE, quote = FALSE)
