#' Generates a matrix of count values which is used to create the ligand occupancy map
#'
#' see ?readCOM() to understand formatting of .xyz file.
#'
#' @param directory Path containing the input .xyz file(s)
#' @param pmat A 4 x 4 Viewing Matrix used to transform the 3D coordinates to 2D. Callend by readCOM() function
#' @param spacing the size of the returend matrix ( 75 * 75 is the default)
#' @return Matrix of count values used for heatmap generation
#'
#' Example
#'
#' pmat <- matrix(rep(0, 16), ncol = 4)
#' diag(pmat) <- 1
#' lom <- lom_matrix(directory = "~/Documents/xyz_files/", pmat = view_mat, spacing = 75)
#' 
#' plot_lom(lom)
#'


lom_matrix <- function(directory, pmat, spacing=75)
{
	XY <- readCOM(directory, pmat)
	#
	Ln    <- length(XY)+1
	x.min <- XY$info[1]; x.max <- XY$info[2]
	y.min <- XY$info[3]; y.max <- XY$info[4]
	#
	iL <- XY$info[5]
	#
	x.seq <- seq(x.min, x.max, abs(x.max-x.min)/spacing)
	y.seq <- seq(y.min, y.max, abs(y.max-y.min)/spacing)
	#
	lom <- matrix(0, ncol=spacing+4, nrow=spacing+4)
	#
	x.mat <- matrix(0, ncol=3, nrow=spacing); colnames(x.mat) <- c("From", "To", "Counter")
	y.mat <- matrix(0, ncol=3, nrow=spacing); colnames(y.mat) <- c("From", "To", "Counter")
	#
	x.mat[,1] <- x.seq[1:(spacing)]
	x.mat[,2] <- x.seq[2:(spacing+1)]
	y.mat[,1] <- y.seq[1:(spacing)]
	y.mat[,2] <- y.seq[2:(spacing+1)]
	#
	for(lLl in 1:iL)
	{
		for(i in 1:length(XY[[lLl]]$x))
		{
			xcount = 0
			ycount = 0
			for(j in 1:spacing)
			{
				if(XY[[lLl]]$x[i] >= x.mat[j,1] && XY[[lLl]]$x[i] < x.mat[j,2])
				{
					xtemp.ind <- j
					xcount = 1
				}
				if(XY[[lLl]]$y[i] >= y.mat[j,1] && XY[[lLl]]$y[i] < y.mat[j,2])
				{
					ytemp.ind <- j
					ycount = 1
				}
				if(xcount == 1 && ycount ==1)
				{
					xtemp.ind <- xtemp.ind + 2
					ytemp.ind <- ytemp.ind + 2
					# +3
					lom[xtemp.ind, ytemp.ind] <- lom[xtemp.ind, ytemp.ind] + 3
					# +2
					lom[xtemp.ind-1, ytemp.ind]   <- lom[xtemp.ind-1, ytemp.ind]   + 2
					lom[xtemp.ind+1, ytemp.ind]   <- lom[xtemp.ind+1, ytemp.ind]   + 2
					lom[xtemp.ind, ytemp.ind-1]   <- lom[xtemp.ind, ytemp.ind-1]   + 2
					lom[xtemp.ind, ytemp.ind+1]   <- lom[xtemp.ind, ytemp.ind+1]   + 2
					lom[xtemp.ind-1, ytemp.ind-1] <- lom[xtemp.ind-1, ytemp.ind-1] + 2
					lom[xtemp.ind+1, ytemp.ind-1] <- lom[xtemp.ind+1, ytemp.ind-1] + 2
					lom[xtemp.ind+1, ytemp.ind+1] <- lom[xtemp.ind+1, ytemp.ind+1] + 2
					lom[xtemp.ind-1, ytemp.ind+1] <- lom[xtemp.ind-1, ytemp.ind+1] + 2
                    # +1
					lom[xtemp.ind-2, ytemp.ind-2] <- lom[xtemp.ind-2, ytemp.ind-2] + 1
					lom[xtemp.ind-1, ytemp.ind-2] <- lom[xtemp.ind-1, ytemp.ind-2] + 1
					lom[xtemp.ind, ytemp.ind-2]   <- lom[xtemp.ind, ytemp.ind-2]   + 1
					lom[xtemp.ind+1, ytemp.ind-2] <- lom[xtemp.ind+1, ytemp.ind-2] + 1
					lom[xtemp.ind+2, ytemp.ind-2] <- lom[xtemp.ind+2, ytemp.ind-2] + 1
					lom[xtemp.ind-2, ytemp.ind-1] <- lom[xtemp.ind-2, ytemp.ind-1] + 1
					lom[xtemp.ind+2, ytemp.ind-1] <- lom[xtemp.ind+2, ytemp.ind-1] + 1
					lom[xtemp.ind-2, ytemp.ind]   <- lom[xtemp.ind-2, ytemp.ind]   + 1
					lom[xtemp.ind+2, ytemp.ind]   <- lom[xtemp.ind+2, ytemp.ind]   + 1
					lom[xtemp.ind-2, ytemp.ind+1] <- lom[xtemp.ind-2, ytemp.ind+1] + 1
					lom[xtemp.ind+2, ytemp.ind+1] <- lom[xtemp.ind+2, ytemp.ind+1] + 1
					lom[xtemp.ind-2, ytemp.ind+2] <- lom[xtemp.ind-2, ytemp.ind+2] + 1
					lom[xtemp.ind-1, ytemp.ind+2] <- lom[xtemp.ind-1, ytemp.ind+2] + 1
					lom[xtemp.ind, ytemp.ind+2]   <- lom[xtemp.ind, ytemp.ind+2]   + 1
					lom[xtemp.ind+1, ytemp.ind+2] <- lom[xtemp.ind+1, ytemp.ind+2] + 1
					lom[xtemp.ind+2, ytemp.ind+2] <- lom[xtemp.ind+2, ytemp.ind+2] + 1
                    print(i)
                    break
				}
			}
		}
	}
	return(lom)
}
