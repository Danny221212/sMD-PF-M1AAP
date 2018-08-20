#' Calculates the point at which the ligand has the strongest ineteraction energy with the protein
#'
#'
#'
#' @param ene_data data frame of combined interaction energy. derived from combine_energies() function.
#' @param files character vector contatining the full paths to the st2 files, identical to the one used in combine_energies() function.
#' @param frames number of frames in the input simulations for energy calculations (or rows containted in the .st2 files).
#' @param lignodes number of requested ligand nodes for the lig-path.
#' @param energy_threshold The energy threshold to included residues in the network.
#' @return data frame contatining the From, TO and Energy columns which an igraph object can be created from
#' @export
#'


energy_relations <- function(ene_data, files, frames = 1250, lignodes = 50, energy_threshold = 100){
    #
    data     <- data.frame()
    resnames <- rownames(ene_data)
    RowStep  <- lignodes/2
    #
    connections <- plyr::ldply(files, .progress = "text", function(path){
                        data     <- data.frame()
                        st2temp <- read_st2(path) 
                        st2temp <- st2temp[[1]] + st2temp[[2]]
                        for(i in 1:(frames/(lignodes/2))){
                            data <- rbind(data, colSums(st2temp[i:(i+lignodes),]))
                        }
                        colnames(data) <- resnames
                        min_col <- apply(data, 2, function(x) which.min(x))
                        return(min_col)
    })
    #
    ligsCons <- round(colMeans(connections), 0)
    data <- data.frame( from = resnames, to = as.vector(as.numeric(ligsCons)), energy = (abs(ene_data[,1])))
    indFORocc <- which(data$energy>1000)
    data <- data[which(data$energy>1000),]
    return(data)
}
