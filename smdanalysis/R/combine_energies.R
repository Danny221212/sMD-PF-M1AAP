#' Reads in all st2 files and combines them into one aggregate data frame.
#'
#' Ensure each st2 file contatins the same number of rows OR the multiple simulations have the same number of frames. 
#'
#' @param files A character vector contatining the full paths to the ST2 files per ligand.
#' @return data frame of per residue interaction energy
#' @export
#'

combine_energies <- function(files){
    ene   <- plyr::ldply(files, .progress = "text", function(path){
            st2temp <- read_st2(path) 
            st2temp <- st2temp[[1]] + st2temp[[2]]
            st2temp <- colSums(st2temp)
            return(st2temp)
    })
    ene <- colMeans(ene)
    return(as.data.frame(ene))
}
