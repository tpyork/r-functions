# ========================================================================
# CREATED: 2015-05-21                                    
# Functions for Illumina 450k beadchip                                     
# ========================================================================

# ------------------------------------------------------------------------
# CONTENTS

# 0.  DMRfind.same()   # Has the restriction that coefs are in the same direction
# 1.  DMRfind()        # No restriction on coef direction
# 2.  DMRfind.v1()     # original function
# 3.  DMRauc()
# 4.  filterCrossReactiveProbes()
# 5.  snpCheckMySamples()


# ------------------------------------------------------------------------
# DMRfindSame() ---------------------------------------------------------
# ------------------------------------------------------------------------

DMRfind.same <- function(dat, k) {
  
  # Find overlapping positions defined by k
  # Has the restriction that coefs are in the same direction
  # data must be ordered by chromosome then position; done within function
  # columns are c('cg','chr','position','statistic','coef')
  # chr is a single digit ('1', '2', 'X')
  # suggest filter cg's by univariate statistic (p,q, 5/95 percentile)
  # dat is data.frame; k is max distance between CpGs
  
  count <- 1
  out <- data.frame(region=NA, chr=NA, start=NA, stop=NA, length=NA, num_cpg=NA, pos_coef=NA, neg_coef=NA)
  clean <- TRUE
  chr.curr <- 1
  
  # sort temp
  temp <- dat[,c('cg','chr','position','statistic','coef')]
  #temp$chr <- as.numeric(temp$chr)
  temp <- temp[order(temp$chr, temp$position),]
  
  
  for (i in 1:dim(temp)[1]) {
    #initialize a new region here
    
    # start a new DMR
    if (clean) {
      new <- data.frame(region=NA, chr=NA, start=NA, stop=NA, length=NA, num_cpg=NA, pos_coef=0, neg_coef=0)    
      cpgs <- 1
      new$region <- count
      new$chr <- temp[i,'chr']
      new$start <- temp[i,'position']
      new$stop <- temp[i,'position']
      new$length <- new$stop - new$start
      new$num_cpg <- 1
      new$pos_coef <- new$pos_coef + ifelse(temp[i,'coef'] >= 0, 1, 0)
      new$neg_coef <- new$neg_coef + ifelse(temp[i,'coef'] < 0, 1, 0)  
      point <- temp[i,'position']
    }
    
    # end if last CpG
    if (i==dim(temp)[1]) {
      out <- rbind(out, new)
      out <- out[-1,]
      break
    }
    
    # end if next CpG on a new chr
    if (temp[i+1,'chr'] != chr.curr) {
      out <- rbind(out, new)    
      clean <- TRUE
      count <- count + 1
      chr.curr <- temp[i+1,'chr']
      next
    }
    
    
    ifelse(is.na(temp[i,'coef']), 0, sign(temp[i,'coef']))
    
    
    # if (within k & same coef direction) then add to current DMR; if not then write record to out
    if ( (point+k >= temp[i+1,'position']) & 
         (ifelse(is.na(temp[i,'coef']), 0, sign(temp[i,'coef']))==ifelse(is.na(temp[i+1,'coef']), 0, sign(temp[i+1,'coef'])) ) ) {
      new$stop <- temp[i+1,'position']
      new$length <- new$stop - new$start
      new$num_cpg <- new$num_cpg + 1
      new$pos_coef <- new$pos_coef + ifelse(temp[i+1,'coef'] >= 0, 1, 0)
      new$neg_coef <- new$neg_coef + ifelse(temp[i+1,'coef'] < 0, 1, 0)  
      point <- new$stop
      clean <- FALSE
      
    } else {
      out <- rbind(out, new)    
      clean <- TRUE
      count <- count + 1
    }
  }
  
  return(out)
}
# ------------------------------------------------------------------------




# ------------------------------------------------------------------------
# DMRfind() --------------------------------------------------------------
# ------------------------------------------------------------------------

DMRfind <- function(dat, k) {

	# Find overlapping position
	# data must be ordered by chromosome then position; done within function
	# columns are c('cg','chr','position','statistic')
	# suggest filter cg's by univariate statistic (p,q, 5/95 percentile)
	# dat is data.frame; k is max distance between CpGs

  count <- 1
  out <- data.frame(region=NA, chr=NA, start=NA, stop=NA, length=NA, num_cpg=NA, pos_coef=NA, neg_coef=NA)
  clean <- TRUE
  chr.curr <- 1

  # sort temp
  temp <- dat[,c('cg','chr','position','statistic','coef')]
  temp <- temp[order(temp$chr, temp$position),]
  #names(temp) <- c('cg','chr','position','statistic','coef')

  for (i in 1:dim(temp)[1]) {
    #initialize a new region here

    # start a new DMR
    if (clean) {
      new <- data.frame(region=NA, chr=NA, start=NA, stop=NA, length=NA, num_cpg=NA, pos_coef=0, neg_coef=0)    
      cpgs <- 1
      new$region <- count
      new$chr <- temp[i,'chr']
      new$start <- temp[i,'position']
      new$stop <- temp[i,'position']
      new$length <- new$stop - new$start
      new$num_cpg <- 1
      new$pos_coef <- new$pos_coef + ifelse(temp[i,'coef'] >= 0, 1, 0)
      new$neg_coef <- new$neg_coef + ifelse(temp[i,'coef'] < 0, 1, 0)  
      point <- temp[i,'position']
    }

    # end if last CpG
    if (i==dim(temp)[1]) {
      out <- rbind(out, new)
      out <- out[-1,]
      break
    }

    # end if next CpG on a new chr
    if (temp[i+1,'chr'] != chr.curr) {
      out <- rbind(out, new)    
      clean <- TRUE
      count <- count + 1
      chr.curr <- temp[i+1,'chr']
      next
    }

    # if within k then add to current DMR; if not then write record to out
    if (point+k >= temp[i+1,'position']) {
      new$stop <- temp[i+1,'position']
      new$length <- new$stop - new$start
      new$num_cpg <- new$num_cpg + 1
      new$pos_coef <- new$pos_coef + ifelse(temp[i+1,'coef'] >= 0, 1, 0)
      new$neg_coef <- new$neg_coef + ifelse(temp[i+1,'coef'] < 0, 1, 0)  
      point <- new$stop
      clean <- FALSE
  
    } else {
      out <- rbind(out, new)    
      clean <- TRUE
      count <- count + 1
    }
  }

  return(out)
}
# ------------------------------------------------------------------------




# ------------------------------------------------------------------------
# DMRfind.v1() -----------------------------------------------------------
# ------------------------------------------------------------------------

DMRfind.v1 <- function(dat, k) {

    # NB: does not include DMR summary of positive or negative coefficients
	# Find overlapping position
	# data must be ordered by chromosome then position
	# columns are c('cg','chr','position','statistic')
	# suggest filter cg's by univariate statistic (p,q, 5/95 percentile)
	# dat is data.frame; k is max distance between CpGs

  count <- 1
  out <- data.frame(region=NA, chr=NA, start=NA, stop=NA, length=NA, num_cpg=NA)
  clean <- TRUE
  chr.curr <- 1

  # sort temp
  temp <- dat[,c('cg','chr','position','statistic')]
  temp <- temp[order(temp$chr, temp$position),]

  for (i in 1:dim(temp)[1]) {
    #initialize a new region here

    # start a new DMR
    if (clean) {
      new <- data.frame(region=NA, chr=NA, start=NA, stop=NA, length=NA, num_cpg=NA)    
      cpgs <- 1
      new$region <- count
      new$chr <- temp[i,'chr']
      new$start <- temp[i,'position']
      new$stop <- temp[i,'position']
      new$length <- new$stop - new$start
      new$num_cpg <- 1
      point <- temp[i,'position']
    }

    # end if last CpG
    if (i==dim(temp)[1]) {
      out <- rbind(out, new)
      out <- out[-1,]
      break
    }

    # end if next CpG on a new chr
    if (temp[i+1,'chr'] != chr.curr) {
      out <- rbind(out, new)    
      clean <- TRUE
      count <- count + 1
      chr.curr <- temp[i+1,'chr']
      next
    }

    # if within k then add to current DMR; if not then write record to out
    if (point+k >= temp[i+1,'position']) {
      new$stop <- temp[i+1,'position']
      new$length <- new$stop - new$start
      new$num_cpg <- new$num_cpg + 1
      point <- new$stop
      clean <- FALSE
  
    } else {
      out <- rbind(out, new)    
      clean <- TRUE
      count <- count + 1
    }
  }

  return(out)
}
# ------------------------------------------------------------------------




# ------------------------------------------------------------------------
# DMRauc() ---------------------------------------------------------------
# ------------------------------------------------------------------------

DMRauc <- function(temp) {

	# Calculate area of multiple trapezoids
	# temp contains only cg's that make up a single DMR
	# Columns of temp are: c('cg','chr','position','statistic')
	# DMRs can be identified from a relaxed threshold of univariate tests

  temp <- temp[order(temp$position),]
  auc.out <- 0

  for (i in dim(temp)[1]:1) {
    h <- temp[i,'position'] - temp[i-1,'position']
    auc <- round(0.5*(temp[i,'statistic'] + temp[i-1,'statistic'])*h,4)
    auc.out <- auc.out + auc
    if (i==2) break
  }

  return(auc.out)
}
# ------------------------------------------------------------------------




# ------------------------------------------------------------------------
# filterCrossReactiveProbes() --------------------------------------------
# ------------------------------------------------------------------------

filterCrossReactiveProbes <- function(gset) {
  #INPUT: RGset or GenomicRatioSet object
  #OUTPUT: same object with the 29233 cross-reactive probes removed; identified in PMID: : 23314698 
  # crp <- read.table("/Users/timothyyork/Documents/Projects/methylation450k/crossReactiveProbes_450k/crossReactiveProbes_450k.txt", sep="\t", header=T)
  # save(crp, file="/Users/timothyyork/Dropbox/R_FXNs/support_files/crossReactiveProbes_450k.Rdata")

  load("/Users/timothyyork/Dropbox/R_FXNs/support-files/crossReactiveProbes_450k.Rdata")

  gset2 <- gset[!is.element(rownames(gset), crp$TargetID), ]
  gset2
}
# ------------------------------------------------------------------------




# ------------------------------------------------------------------------
# snpCheckMySamples-------------------------------------------------------
# ------------------------------------------------------------------------

snpCheckMySamples <- function(RGset_object, twins = FALSE, repeated_measure = FALSE, common_ID = NULL, 
                              individual_ID, highMatchThreshold = .95, lowMatchThreshold = .65){
  
  # Pass in an RGset and column name that has (at least) individual identifiers in the pData.
  # If this sample includes repeated measures or twins, common_ID can be used for the common or family identifier.
  
  # The purpose of this function is to ensure that cryptic duplicate samples do not exist in the RGset. 
  # If the sample includes multiple measurements from the same individual (or twins), use the "common_ID" argument
  # to specify the family identifier vs. individual identifier (individual_ID).
  # Additionally, this function can check the zygosity of twins.
  
  
  ind_index <- which(names(pData(RGset_object)) == individual_ID)
  fam_index <- which(names(pData(RGset_object)) == common_ID)
  
  
  if(length(ind_index) == 0){
    
    return("There appears to be a problem with the individual_ID. Please check your spelling and your object.")
    
  }
  
  
  # Repeated samples and/or twins -----
  if (twins == TRUE | repeated_measure == TRUE) {
    
    
    if(length(fam_index) == 0){
      
      return("There appears to be a problem with the common_ID. Without this information, zygosity cannot be checked. Please check your spelling and your object.")
      
    }
    
    
    mySnpBeta <- getSnpBeta(RGset_object)
    tempNames <- paste0(pData(RGset_object)[ ,ind_index], "/" , pData(RGset_object)[ ,fam_index])
    colnames(mySnpBeta) <- tempNames # Adds sample & family names
    snp.cors <- cor(mySnpBeta)                         # make correlation matrix
    snp.cors[upper.tri(snp.cors, diag = TRUE)] <- NA   # remove duplicates
    
    
    cor.df <- as.data.frame.table(snp.cors, stringsAsFactors = FALSE)
    names(cor.df) <- c("ID1", "ID2", "Corr")
    cor.df <- subset(cor.df, !is.na(Corr))
    
    
    
    
    # Create Histogram of Relatedness---------
    hist(cor.df$Corr, breaks = "Scott", col = "dark grey", main = "Histogram of Distribution of Correlated SNP Beta Values")
    
    
    
    
    # Separate Individual and Family IDs---------
    cor.df$fam_ID1 <-  sapply(strsplit(cor.df$ID1, "/"), `[`, 2) # family
    cor.df$fam_ID2 <-  sapply(strsplit(cor.df$ID2, "/"), `[`, 2)
    
    cor.df$ID1 <- sapply(strsplit(cor.df$ID1, "/"), `[`, 1) # individual
    cor.df$ID2 <- sapply(strsplit(cor.df$ID2, "/"), `[`, 1) # individual
    
    
    
    
    # Identify samples that do not match closely enough---------
    sameID <- cor.df[which(cor.df$fam_ID1 == cor.df$fam_ID2),]
    problem1 <- sameID[which(sameID$Corr < highMatchThreshold), ]
    
    
    
    # Identify samples that match too closely---------
    problem2 <- cor.df[which(cor.df$fam_ID1 != cor.df$fam_ID2 & cor.df$Corr > lowMatchThreshold), ]
    
    alist <- list(corrs = cor.df, lowMatch = problem1, highMatch = problem2)
    
    return(alist)
  }
  
  
  
  # Individual Samples only ----
  
  else{
    
    mySnpBeta <- getSnpBeta(RGset_object)
    colnames(mySnpBeta) <- pData(RGset_object)[ ,ind_index]  # Adds sample names
    snp.cors <- cor(mySnpBeta)                         # make correlation matrix
    snp.cors[upper.tri(snp.cors, diag = TRUE)] <- NA   # remove duplicates
    
    
    cor.df <- as.data.frame.table(snp.cors, stringsAsFactors = FALSE)
    names(cor.df) <- c("ID1", "ID2", "Corr")
    cor.df <- subset(cor.df, !is.na(Corr))
    
    
    # Create Histogram of Relatedness---------
    hist(cor.df$Corr, breaks = "Scott", col = "dark grey", main = "Histogram of Distribution of Correlated SNP Beta Values")
    
    
    
    
    
    # Identify samples that match too closely
    problem2 <- cor.df[which(cor.df$ID1 != cor.df$ID2 & cor.df$Corr > lowMatchThreshold), ]
    alist <- list(corrs = cor.df, highMatch = problem2)
    
    return(alist)
  }
  
}








