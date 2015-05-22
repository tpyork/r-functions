# ========================================================================
# CREATED: 05.21.15                                     
# Find DMR                                     
# ========================================================================



# SET UP WORKSPACE AND LOAD PACKAGES -------------------------------------

setwd("/Users/tpyork/Documents/Projects/CaB/")
load("CaB_Pre_Post_Me.Rdata")





# SIMULATE DATA ----------------------------------------------------------

temp <- data.frame(cg= paste('cg',1:5,sep=''), chr= rep(1,5), position= c(2000,5000,5500,5525,8000), statistic= c(8,3.6, 5, 5.5, 3))



# FIND OVERLAPPING POSITIONS ---------------------------------------------

k <- 300
count <- 1
out <- data.frame(region=NA, chr=NA, start=NA, stop=NA, length=NA, num_cpg=NA)
clean <- TRUE

for (i in 1:dim(temp)[1]) {
  #initialize a new region here
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

  if (i==dim(temp)[1]) {
    out <- rbind(out, new)
    out <- out[-1,]
    break
  }

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



















