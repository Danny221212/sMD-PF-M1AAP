#' Reads in a PDB file and converts atomic coordinates to residue centroids XYZ.
#'
#'
#'
#' @param pdb protein data bank file
#' @return data frame of per residue centroids containing Residue Number, Residue Name, x, y and z coordinates. 
#' @export
#'


pdb_res_to_centroid <- function(pdb = "input.pdb"){
    #
    protein 	<- read.fwf(pdb, widths=c(6,6,4,1,4,1,4,4,8,8,8,6,12,4,2,2), as.is = T, stringsAsFactors = FALSE, header=FALSE)
    protein 	<- protein[-which(! "ATOM  " == as.vector(as.character(protein[,1]))),]
    protein     <- apply(protein, 2, function(x) gsub('\\s+', '',x))
    protein     <- as.data.frame(protein)
    protein$V2  <- as.numeric(as.character(protein$V2));  protein$V7   <- as.numeric(as.character(protein$V7))
    protein$V8  <- as.numeric(as.character(protein$V8));  protein$V9   <- as.numeric(as.character(protein$V9))
    protein$V10 <- as.numeric(as.character(protein$V10)); protein$V11  <- as.numeric(as.character(protein$V11))
    protein$V12 <- as.numeric(as.character(protein$V12));  protein$V13 <- as.numeric(as.character(protein$V13))
    #
    colnames(protein) <- c("Type","AtonNumber","AtomName","AltLoc","ResName", "ChainName", "ResNumber", "Insert", "x", "y", "z", "Occupancy", "BFactor", "SEGid", "Element", "Charge")
    #
    protein 	<- protein[-length(protein[,1]),]
    resideTOT 	<- unique(protein$ResNumber)
    #
    xMean <-aggregate(as.numeric(as.vector(protein$x)), by=list(ResNumber=as.numeric(protein$ResNumber)), FUN=mean)
    yMean <-aggregate(as.numeric(as.vector(protein$y)), by=list(ResNumber=as.numeric(protein$ResNumber)), FUN=mean)
    zMean <-aggregate(as.numeric(as.vector(protein$z)), by=list(ResNumber=as.numeric(protein$ResNumber)), FUN=mean)
    #
    protein$ResName <- as.vector(as.character(protein$ResName))
    resname <- vector()
    for(i in 1:length(resideTOT)){
        resname[i] <- protein[which(protein$ResNumber==resideTOT[i])[1], "ResName"]
    }
    REScentroid <- cbind(resname, xMean, yMean[,2], zMean[,2])
    colnames(REScentroid) <- c("resname", "resnumber", "x","y","z")
    return(REScentroid)
}

