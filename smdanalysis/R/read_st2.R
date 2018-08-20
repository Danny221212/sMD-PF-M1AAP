#' Reads in st2 files in a given directory
#'
#'
#'
#' @param st2dir single directory contatining all st2 files
#' @param as.sum Default is FALSE which seperates VdW and Coloumb Energies
#' @return List of interaction energies.
#' @export
#'


read_st2 <- function(st2dir, as.sum=FALSE)
{
 	cwd   <- getwd(); setwd(st2dir)
	files <- list.files() #list.files("st2dir") 
	ind   <- grep("out", files)
	files <- files[ind]
	files <- mixedsort(files) #list.files("st2dir")
	resNAMES <- vector()
	counter = 1
	for(i in files)
	{
		info  <- file.info(i)
		if(info$size == 0){
			cat("WARNING: File",i,"is Empty\n")
			next
		}
		if(i == files[1])
		{
			data   <- read.delim(i, skip="6", sep=" ")
			vdw    <- as.numeric(data[,12])
			columb <- gsub("\\[", "", data[,11]); columb <- as.numeric(columb)
			indC   <- which(is.na(columb)==TRUE); columb <- columb[-indC]
			indV   <- which(is.na(vdw)==TRUE); vdw <- vdw[-indV]
			vdwMAT <- matrix(vdw)	
			colMAT <- matrix(columb)
		} else {
			data   <- read.delim(i, skip="6", sep=" ")
			vdw    <- as.numeric(data[,12])
			columb <- gsub("\\[", "", data[,11]); columb <- as.numeric(columb)
			indC   <- which(is.na(columb)==TRUE); columb <- columb[-indC]
			indV   <- which(is.na(vdw)==TRUE); vdw <- vdw[-indV]
			colMAT <- cbind(colMAT, columb)
			vdwMAT <- cbind(vdwMAT, vdw)
		}
		resNAMES[counter] <- i
		counter <- counter+1
	}
	resNAMES <- gsub("L_", "", resNAMES); resNAMES <- gsub("S_", "", resNAMES); resNAMES <- gsub("_out", "", resNAMES); resNAMES <- gsub("RE", "", resNAMES)
	colnames(colMAT) <- resNAMES; colnames(vdwMAT) <- resNAMES
	setwd(cwd)
	if(as.sum==FALSE)
	{
		l1 <- list(colMAT, vdwMAT); names(l1) <- c("columb", "vdw")
		return(l1)
	} else {
		colMAT <- apply(colMAT, 2, sum);vdwMAT <- apply(vdwMAT, 2, sum)
		l1 <- matrix(colMAT, nrow=1); l1 <- rbind(colMAT, vdwMAT)
		rownames(l1)  <- c("columb", "vdw")
		return(l1)
	}
}
