#' Loads multiple or a single Center Of Mass (COM) or Centroid .xyz file(s) and 
#' Converts to 2D x and y coordinates based on a viewing matrix
#'
#' NB*** Ensure the extension of the files is .xyz, contains no header and is column ordered X, Y, Z.
#'
#' @param directory Path containing the input .xyz file(s)
#' @param pmat A 4 x 4 Viewing Matrix used to transform the 3D coordinates to 2D 
#' @param seperator how the columns are seperated. Example: " " for space, "," for csv files, "\t" for tab-delim files.
#' @return list with transformed coordinates and min/max + length of files loaded info
#' @export
#'
#'
#' Example
#'
#' pmat <- matrix(rep(0, 16), ncol = 4)
#' diag(pmat) <- 1
#' lom <- lom_matrix(directory = "~/Documents/xyz_files/", pmat = view_mat, spacing = 75)
#' 
#' plot_lom(lom)
#'

readCOM <- function(directory, pmat, seperator = " ", extension = ".xyz")
{
	files <- list.files(directory, pattern = extension)
	inputLength <- length(files)
	if(inputLength == 0){
		stop("No *.xyz files in directory.")
	}
	output <- list()
	min_x  <- 0
	min_y  <- 0
	max_x  <- 0
	max_y  <- 0
	for(i in 1:inputLength){
		dataTemp <- read.delim(files[i], sep = seperator, header = FALSE)
		colnames(dataTemp) <- c("x", "y", "z")
		dataTemp <- trans3D(dataTemp$x, dataTemp$y, dataTemp$z, pmat = pmat) 
		output[[i]] <- dataTemp
    }
	if(min_x == 0){
		min_x <- min(dataTemp$x)
		max_x <- max(dataTemp$x)
		min_y <- min(dataTemp$y)
		max_y <- max(dataTemp$y)
	}
	# MIN #
	if(min_x > min(dataTemp$x)){
		min_x <- min(dataTemp$x)
	}
	if(min_y > min(dataTemp$y)){
		min_y <- min(dataTemp$y)
	}
	# MAX #
	if(max_x < max(dataTemp$x)){
		max_x <- max(dataTemp$x)
	}
	if(max_y < max(dataTemp$y)){
		max_y <- max(dataTemp$y)
	}
	output$info <- c(min_x, max_x, min_y, max_y, inputLength)
	return(output)
}
