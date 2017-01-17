#########################################################################
#targets <- read.450k.sheet(baseDir)
#RGset <- read.450k.exp(base = baseDir, targets = targets,extended=T)
#Beta = getBeta(preprocessRaw2(RGset))
#SNP.probes = Beta[grep("rs",rownames(Beta)),]


################################################################################
getManifestInfo2 <- function(object, type = c("nLoci", "locusNames")) {
	type <- match.arg(type)
	switch(type,
			"nLoci" = {
				nrow(getProbeInfo(object, type = "I")) +
						nrow(getProbeInfo(object, type = "SnpI")) +
						nrow(getProbeInfo(object, type = "II"))+ 
						nrow(getProbeInfo(object, type = "SnpII"))
			},
			"locusNames" = {
				c(getProbeInfo(object, type = "I")$Name,
						getProbeInfo(object, type = "SnpI")$Name,
							getProbeInfo(object, type = "II")$Name,
							getProbeInfo(object, type = "SnpII")$Name)
			})
}




preprocessRaw2 <- function(rgSet) {
	if(!is(rgSet, "RGChannelSet"))
		stop("'rgSet' needs to be a 'RGChannelSet'")
	locusNames <- getManifestInfo2(rgSet, type="locusNames")
	M <- matrix(NA_real_, ncol = ncol(rgSet), nrow = length(locusNames),
			dimnames = list(locusNames, sampleNames(rgSet)))
	U <- M
	TypeII.Name <- getProbeInfo(rgSet, type = "II")$Name
	TypeSnpII.Name <- getProbeInfo(rgSet, type = "SnpII")$Name
	
	
	M[TypeII.Name,] <- getGreen(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA,]
	U[TypeII.Name,] <- getRed(rgSet)[getProbeInfo(rgSet, type = "II")$AddressA,]
	
	M[TypeSnpII.Name,] <- getGreen(rgSet)[getProbeInfo(rgSet, type = "SnpII")$AddressA,]
	U[TypeSnpII.Name,] <- getRed(rgSet)[getProbeInfo(rgSet, type = "SnpII")$AddressA,]
	
	
	
	TypeI.Red <- getProbeInfo(rgSet, type = "I-Red")
	TypeI.Green <- getProbeInfo(rgSet, type = "I-Green")
	
	SNPI <- getProbeInfo(rgSet, type = "SnpI")
	
	
	M[TypeI.Red$Name,] <- getRed(rgSet)[TypeI.Red$AddressB,]
	M[TypeI.Green$Name,] <- getGreen(rgSet)[TypeI.Green$AddressB,]
	U[TypeI.Red$Name,] <- getRed(rgSet)[TypeI.Red$AddressA,]
	U[TypeI.Green$Name,] <- getGreen(rgSet)[TypeI.Green$AddressA,]
	
	M[SNPI$Name,] <- getRed(rgSet)[SNPI$AddressB,]
	M[SNPI$Name,] <- getGreen(rgSet)[SNPI$AddressB,]
	U[SNPI$Name,] <- getRed(rgSet)[SNPI$AddressA,]
	U[SNPI$Name,] <- getGreen(rgSet)[SNPI$AddressA,]
	
	
	out <- new("MethylSet",Meth = M, Unmeth = U,
			phenoData = phenoData(rgSet),
			annotation = annotation(rgSet))
	
	out@preprocessMethod <- c(rg.norm = "Raw (no normalization or bg correction)",
			minfi = as.character(packageVersion("minfi")),
			manifest = as.character(packageVersion("IlluminaHumanMethylation450kmanifest")))
	#packageVersion(getManifest(rgSet)))
	out
}
################################################################################