# ========================================================================
# CREATED: 05.21.15                                     
# Find DMR                                     
# ========================================================================



# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------

setwd("/Users/tpyork/Documents/Projects/CaB/")
load("CaB_Pre_Post_Me.Rdata")





# SIMULATE DATA ----------------------------------------------------------


temp <- data.frame(cg= paste('cg',1:7,sep=''), chr= c(rep(1,7)), position= c(2000,5000,5500,5525,8000,8100,8500), statistic= c(8,3.6, 5, 5.5, 3, 4,4.5))

temp <- data.frame(cg= paste('cg',1:12,sep=''), chr= c(rep(1,7), rep(2,4), rep('X',1)), position= c(2000,5000,5500,5525,8000,8100,8500, 2000,5000,5500,9000,9100), statistic= c(8,3.6, 5, 5.5, 3, 4,4.5, 3, 4, 2.5, 7, 5))

temp <- data.frame(cg= paste('cg',1:4,sep=''), chr= c(rep(1,2), rep(2,2)), position= c(2000,5000,3000,3500), statistic= c(3, 4, 2.5, 7))



# FIND OVERLAPPING POSITIONS ---------------------------------------------
# data must be ordered by chromosome then position
# columns are c('cg','chr','position','statistic')
# suggest filter cg's by univariate statistic (p,q, 5/95 percentile)


DMRfind <- function(dat, k) {

  #k <- 1000
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

DMRfind(temp, 1000)

















