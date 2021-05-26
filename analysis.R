library(ShortRead)
library(Biostrings)
library(dplyr)
library(ggplot2)

# Read in sample metadata.
samples <- readr::read_delim('sampleConfig.inward.allSamples.tsv', '\t')


# Configuration.
R2.analysis.width <- 100
R2.recognitionSeqWidth <- 32
cutadapt <- '/usr/bin/cutadapt3'
mafft <- '/home/everett/ext/mafft/bin/mafft'
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
suppressWarnings(file.remove(c('R2.ff', 'R2.ff.out', 'R2.log')))


# Quality trim R2 reads.
R2 <- trimTailw(R2, 2, '?', 5)
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
      
       
       # Identify 5' and 3' reads by selecting reads which match the first 50 NT of the two options.
       # Allow 1 mismatch per 10 NTs.
       reads.5p <- reads[vcountPattern(p5.leader, subseq(reads@sread, 1, R2.recognitionSeqWidth), max.mismatch = round(R2.recognitionSeqWidth/10)) > 0]
       reads.5p.ids <- reads.5p@id
       reads.5p <- subseq(reads.5p@sread, 1, R2.analysis.width)
       names(reads.5p) <- reads.5p.ids
       
       reads.3p <- reads[vcountPattern(p3.leader, subseq(reads@sread, 1, R2.recognitionSeqWidth), max.mismatch = round(R2.recognitionSeqWidth/10)) > 0]
       reads.3p.ids <- reads.3p@id
       reads.3p <- subseq(reads.3p@sread, 1, R2.analysis.width)
       names(reads.3p) <- reads.3p.ids
       
       tibble(subject = x$subject, sample = x$sample, replicate = x$replicate,
              reads_5p = list(unique(reads.5p)), reads_3p = list(unique(reads.3p)))
}))


# Collapse replicate reads.
r2 <- bind_rows(lapply(split(r, r$subject), function(x){
       
        p5 <- unique(Reduce('append', x$reads_5p))
        p3 <- unique(Reduce('append', x$reads_3p))
  
        if(length(p5) == 0 | length(p3) == 0) return(data.frame())
        
        p5.d <- stringdist::stringdist(substr(p5.expected.200, 1, R2.analysis.width), p5)
        p3.d <- stringdist::stringdist(substr(p3.expected.200, 1, R2.analysis.width), p3)
        
        # Debug
        #-----------------------------
        p5.expected <- DNAStringSet(substr(p5.expected.200, 1, R2.analysis.width)); names(p5.expected) <- 'Expected'
        p3.expected <- DNAStringSet(substr(p3.expected.200, 1, R2.analysis.width)); names(p3.expected) <- 'Expected'
        
        invisible(lapply(list(c(30, 40), c(41, 50), c(51, 60)), function(d){
          writeFasta(unlist(DNAStringSetList(p5.expected, p5[p5.d >= d[1] & p5.d <= d[2]])), file = paste0('tmp/', x$subject[1], '.5p.dist_', d[1], '-', d[2], '.fasta'))
          writeFasta(unlist(DNAStringSetList(p3.expected, p3[p3.d >= d[1] & p3.d <= d[2]])), file = paste0('tmp/', x$subject[1], '.3p.dist_', d[1], '-', d[2], '.fasta'))
          system(paste0(mafft, ' --clustalout --maxiterate 1000 --localpair tmp/', x$subject[1], '.5p.dist_', d[1], '-', d[2], '.fasta > tmp/', x$subject[1], '.5p.dist_', d[1], '-', d[2], '.aln'))
          system(paste0(mafft, ' --clustalout --maxiterate 1000 --localpair tmp/', x$subject[1], '.3p.dist_', d[1], '-', d[2], '.fasta > tmp/', x$subject[1], '.3p.dist_', d[1], '-', d[2], '.aln'))
        }))
        #------------------------------
        
        a <- data.frame(table(p5.d))
        names(a) <- c('distance', 'count')
        a$source <- 'p5'
        
        b <- data.frame(table(p3.d))
        names(b) <- c('distance', 'count')
        b$source <- 'p3'
        
        a <- bind_rows(a, b)
        a$distance <- as.integer(a$distance)
        a$subject <- x$subject[1]
        a
}))

invisible(file.remove(list.files('tmp', pattern = '*.fasta', full.names = TRUE)))

ggplot(r2, aes(distance, count)) + 
  theme_bw() +
  geom_col() +
  facet_wrap(subject~source, scales = 'free_y', ncol = 2) +
  labs(x = 'Edit distance', y = paste0('Unique ', R2.analysis.width, ' NT reads'))


ggplot(r2, aes(distance, count)) + 
  theme_bw() +
  geom_col() +
  facet_grid(subject~., scales = 'free_y') +
  labs(x = 'Edit distance', y = paste0('Unique ', R2.analysis.width, ' NT reads')) 
