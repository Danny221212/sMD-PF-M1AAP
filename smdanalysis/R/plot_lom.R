#' Plots the ligand occupancy map.
#'
#' The map uses a 6 gradient colour scale with the following colours:
#' (low Count) "darkblue", "cyan", "green", "yellow", "orange", "red" (High Count)
#'
#'
#' @param lom the ligand occupancy map matrix derived using the lom_matrix() function
#' @return NULL; Plots lom as a heatmap
#' @export
#'
#' Example
#'
#' pmat <- matrix(rep(0, 16), ncol = 4)
#' diag(pmat) <- 1
#' lom <- lom_matrix(directory = "~/Documents/xyz_files/", pmat = view_mat, spacing = 75)
#' 
#' plot_lom(lom)
#'

plot_lom <- function(lom){
    if(class(lom) != "matrix"){
        stop("lom must be of class matrix.")
    }
    colfunc <- colorRampPalette(c("darkblue", "cyan", "green", "yellow", "orange", "red"))
    image(lom, col = colfunc(256), axes=FALSE,
        xlab = "", ylab = "",
        main = "Ligand Occupancy Map")
    return(NULL)
}