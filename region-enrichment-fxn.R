# ========================================================================
# CREATED: 2018-08-28                                    
# Functions for Illumina 450k beadchip                                     
# ========================================================================

# ------------------------------------------------------------------------
# CONTENTS


#  regionEnrich()
#  grBootSample()



# ------------------------------------------------------------------------
# DMRenrichPercent() -----------------------------------------------------
# ------------------------------------------------------------------------

regionEnrich <- function(sig, ann, bkg, k, CI) {
  ### *Calculates significance based on percent overlap; not by OR as initial version (DMRenrich)*
  #sig= GRanges object of significant regions of width >= 1
  #ann= GRanges object of annotations to test for enrichment
  #bkg= Background GRanges object for sig and ann
  #k= the number of permutations to estimate significance of enrichment in sig compared to bkg
  #If CI=TRUE then a bootstrap (i= 1000) to estimate 95% CI of sig-ann overlap
  #Depends on grBootSample(); GenomicRanges; doParallel; foreach
  

  #Calculate observed percentage of overlap
  Per.obs <- suppressWarnings(sum(width(intersect(ann, sig, ignore.strand=T)))) / sum(width(sig))
  Per.bkg <- suppressWarnings(sum(width(intersect(ann, bkg, ignore.strand=T)))) / sum(width(bkg))

  #Estimate 95% CIs for proportion overlap using 1000 bootstrap samples
  if(CI) {
      s <- foreach(i=1:1000, .combine= c, .packages="GenomicRanges") %dopar% {
  
      n <- length(sig)
      sig3 <- sig[sample(n, replace= TRUE)]

      #Calculate bootstrap sample percentage of overlap
 
      a <- suppressWarnings(findOverlaps(sig3, ann))
      Per.boot <- sum(width(suppressWarnings(pintersect(sig3[queryHits(a)], ann[subjectHits(a)])))) / sum(width(sig3))
      Per.boot

    }

    Per.boot <- sort(s)
    Per.lower <- Per.boot[0.025 * length(Per.boot)]
    Per.upper <- Per.boot[0.975 * length(Per.boot)]
  } else {
    Per.boot <- NA
    Per.lower <- NA
    Per.upper <- NA
  }

  #Calculate empirical odds ratio distribution
  bkg <- sort(bkg)

  r <- foreach(i=1:k, .combine= c, .packages="GenomicRanges") %dopar% {
    
    #create null samples of GRanges from bkg of similar structure as sig
    sig1 <- as.data.frame(sig)
    sig2 <- data.frame()
    for (i in 1:dim(sig1)[1]) {
      sig2 <- rbind(sig2, grBootSample(sig1[i,], bkg))
    }  

    sig2 <- makeGRangesFromDataFrame(sig2)
    b <- suppressWarnings(findOverlaps(sig2, ann))
    Per.emp <- sum(width(suppressWarnings(pintersect(sig2[queryHits(b)], ann[subjectHits(b)])))) / sum(width(sig2))
    Per.emp

  }

  pval.enrich <- (sum(r > Per.obs) + 1) / (k + 1)
  pval.deplete <- (sum(r < Per.obs) + 1) / (k + 1)
  mean.emp <- mean(r, na.rm=T)
  fold <- Per.obs / mean.emp


  return(list(Per.bkg= Per.bkg, Per.obs= Per.obs, Per.lower= Per.lower, Per.upper= Per.upper, 
              pval.enrich= pval.enrich, pval.deplete= pval.deplete, 
              mean.emp= mean.emp, fold= fold, Per.emp= r, Per.boot= Per.boot))

}




# ------------------------------------------------------------------------
# grBootSample() ---------------------------------------------------------
# ------------------------------------------------------------------------

grBootSample <- function(SIG, BKG) {
  #SIG= a row from a data.frame  
  #BKG= the background GRanges object
  #Returns a random sampling from BKG (a GRanges object) of the same structure as sig1
  
  sig.row <- makeGRangesFromDataFrame(SIG)
  w <- width(sig.row)  
  
  bkg1 <- BKG[width(BKG) >= w]  
  bkg2 <- bkg1[sample(1:length(bkg1), 1)]
  
  starts <- (start(bkg2)):(end(bkg2)-w+1)
  start <- ifelse(length(starts)==1, starts, sample(starts,1))
  as.data.frame(GRanges(seqnames= seqnames(bkg2), 
                        IRanges(start= start, end= start+w-1), 
                        strand= '*'))
}




# EOF --------------------------------------------------------------------
