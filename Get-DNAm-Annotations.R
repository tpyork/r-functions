# ============================================================
# created: 2016-08-26
# Retrieve uselful DNAm annotations and load into working memory
# ============================================================



# SETUP ======================================================

#source("https://bioconductor.org/biocLite.R")

library(GenomicRanges)
library(AnnotationHub)


# GET ANNOTATION DATA FOR CpGi -------------------------------

ah <- AnnotationHub()
ah <- subset(ah, species=="Homo sapiens")
#unique(ah$dataprovider)


# CpG Islands
#query(ah, "CpG")
#mcols(query(ah, "CpG"))
# AH5086; CpG Island UCSC HG19
#ah2 <- display(ah)
# Select CpG Island track to use
CpGi <- ah[["AH5086"]]

#qplot(width(CpGi[width(CpGi) < 5000]))
#mean(width(CpGi[width(CpGi) < 5000]))


# Define North Shore region; 2kb 5' flanking CpGi
temp <- CpGi
Nshore <- GRanges(seqnames= seqnames(temp),
                  ranges= IRanges(start= start(temp)-2001, end= start(temp)-1),
                  strand= '*')


# Define South Shore region; 2kb 3' flanking CpGi
temp <- CpGi
Sshore <- GRanges(seqnames= seqnames(temp),
                  ranges= IRanges(start= end(temp)+1, end= end(temp)+2001),
                  strand= '*')


# Define North Shelf region; 2kb 5' flanking Nshore
Nshelf <- GRanges(seqnames= seqnames(temp),
                  ranges= IRanges(start= start(Nshore)-2001, end= start(Nshore)-1),
                  strand= '*')


# Define South Shelf region; 2kb 3' flanking Sshore
Sshelf <- GRanges(seqnames= seqnames(temp),
                  ranges= IRanges(start= end(Sshore)+1, end= end(Sshore)+2001),
                  strand= '*')




# Get ChromHMM Annotation ------------------------------------

#mcols(query(ah, "ChromHMM"))$genome
#ah2 <- display(ah)

ChromHMM <- ah[["AH46969"]]



# Get Gene Symbols and Gene Granges

data(genesymbol, package= "biovizBase")




# Save version for reproducibility ---------------------------

#save(CpGi,Nshore,Nshelf,Sshore,Sshelf,ChromHMM, file= "/Users/timothyyork/Documents/Projects/CaB/DNAm_CaB/data/AnnotationHub-dataUsed.Rdata")




rm(temp, ah)



# EOF --------------------------------------------------------
