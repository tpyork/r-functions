# ========================================================================
# CREATED: 2015-06-22                                    
# General Functions
# ========================================================================

# ------------------------------------------------------------------------
# CONTENTS

# 1.  lsos()       #reports summary of object size
# 2.  twindat()    #put data from twins on same row



# ------------------------------------------------------------------------
# lsos() -----------------------------------------------------------------
# ------------------------------------------------------------------------

.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
    napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
    names <- ls(pos = pos, pattern = pattern)
    obj.class <- napply(names, function(x) as.character(class(x))[1])
    obj.mode <- napply(names, mode)
    obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
    obj.size <- napply(names, object.size)
    obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
    vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
    obj.dim[vec, 1] <- napply(names, length)[vec]
    out <- data.frame(obj.type, obj.size, obj.dim)
    names(out) <- c("Type", "Size", "Rows", "Columns")
    if (!missing(order.by))
        out <- out[order(out[[order.by]], decreasing=decreasing), ]
    if (head)
        out <- head(out, n)
    out
}
lsos <- function(..., n=100) {
    .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}
# ------------------------------------------------------------------------




# ------------------------------------------------------------------------
# twindat() --------------------------------------------------------------
# ------------------------------------------------------------------------

#FUNCTION TO PUT TWIN PAIRS IN SAME RECORD
#INPUT:  data.frame where each data for each individual is listed in separate record
#OUTPUT: data.frame where twin pairs are listed in a single record
twindat <- function(dat, famid, twinid, zygosity) {
  datA <- dat[dat[,twinid]==min(dat[,twinid]),]    #twin1
  datB <- dat[dat[,twinid]==max(dat[,twinid]),]    #twin2 
  DAT <- merge(datA, datB, by=famid, all.x=TRUE, all.y=TRUE, suffixes=c("_T1","_T2"))
  #DAT[,paste(twinid,"_T2",sep="")] <- NULL
  #DAT[,paste(zygosity,"_T2",sep="")] <- NULL  #
  #DAT[,twinid] <- ifelse(is.na(DAT[,paste(twinid,"_T1",sep="")]), DAT[,paste(twinid,"_T2",sep="")], DAT[,paste(twinid,"_T1",sep="")])
  #DAT[,paste(twinid,"_T1",sep="")] <- NULL
  DAT[,zygosity] <- ifelse(is.na(DAT[,paste(zygosity,"_T1",sep="")]),DAT[,paste(zygosity,"_T2",sep="")],DAT[,paste(zygosity,"_T1",sep="")])
  #DAT[,paste(zygosity,"_T1",sep="")] <- NULL  #
  #DAT$twinid <- NULL
  return(DAT)
}
# ------------------------------------------------------------------------




