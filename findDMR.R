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

k <- 1000
count <- 0
out <- data.frame(region=NA, chr=NA, start=NA, length=NA, num_cpg=NA, area=NA)

for (i in 1:dim(temp)[1]) {
  count <- count + 1
  if (temp[i,'position']+k >= temp[i+1,'position']) {
    count <- count + 1
    out[count,'region'] <- count
    out[count,'chr'] <- temp[i,'chr']
    out[count,
    num_cpgs <- 
  
  } else {
    out[count,'region'] <- count
    out[count,'chr'] <- temp[i,'chr']
    out[count,
    


  }


}



















