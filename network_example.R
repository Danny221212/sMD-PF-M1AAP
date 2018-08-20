# Create iGraph Function
mycircle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.frame.width <- params("vertex", "frame.width")
  if (length(vertex.frame.width) != 1 && !is.null(v)) {
    vertex.frame.width <- vertex.frame.width[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.frame.color,
         vertex.size, vertex.frame.width,
         FUN=function(x, y, bg, fg, size, lwd) {
           symbols(x=x, y=y, bg=bg, fg=fg, lwd=lwd,
                   circles=size, add=TRUE, inches=FALSE)
         })
}

add.vertex.shape("fcircle", clip=igraph.shape.noclip,
		plot=mycircle, parameters=list(vertex.frame.color=1,
                                  vertex.frame.width=4))

######### MAIN BODY #############

setwd("F:/ENERGY_CALCULATIONS/Energy/")

m1.info <- read.csv("F:/ENERGY_CALCULATIONS/ResNAMES.csv", header = TRUE)
ASres   <- c("317", "319", "320","459" , "460", "461", "462", "463", "492", "493", "497", "518", "575", "580", "1034")

output <- list.dirs(path = ".", full.names = TRUE, recursive = TRUE); ST2dirs <- grep("/ST2FILES", output)
output <- output[ST2dirs]


ARG.long.dirs     <- output[grep("ARG",output)]; ARG.long.dirs     		<- ARG.long.dirs[grep("LONG",ARG.long.dirs)]
ARG.short.dirs    <- output[grep("ARG",output)]; ARG.short.dirs     	<- ARG.short.dirs[grep("SHORT",ARG.short.dirs)]
ArgAla.long.dirs  <- output[grep("ArgAla",output)]; ArgAla.long.dirs  	<- ArgAla.long.dirs[grep("LONG",ArgAla.long.dirs)]
ArgAla.short.dirs <- output[grep("ArgAla",output)]; ArgAla.short.dirs  	<- ArgAla.short.dirs[grep("SHORT",ArgAla.short.dirs)]
MetPhe.long.dirs  <- output[grep("MetPhe",output)]; MetPhe.long.dirs  	<- MetPhe.long.dirs[grep("LONG",MetPhe.long.dirs)]
MetPhe.short.dirs <- output[grep("MetPhe",output)]; MetPhe.short.dirs  	<- MetPhe.short.dirs[grep("SHORT",MetPhe.short.dirs)]
BES.long.dirs 	  <- output[grep("BES",output)]; BES.long.dirs     		<- BES.long.dirs[grep("LONG",BES.long.dirs)]
BES.short.dirs 	  <- output[grep("BES",output)]; BES.short.dirs     	<- BES.short.dirs[grep("SHORT",BES.short.dirs)]
R5X.long.dirs 	  <- output[grep("R5X",output)]; R5X.long.dirs     		<- R5X.long.dirs[grep("LONG",R5X.long.dirs)]
R5X.short.dirs 	  <- output[grep("R5X",output)]; R5X.short.dirs     	<- R5X.short.dirs[grep("SHORT",R5X.short.dirs)]


BES.long <- list(); BES.short <- list(); BES.cl = 1; BES.cs = 1
R5X.long <- list(); R5X.short <- list(); R5X.cl = 1; R5X.cs = 1
ARG.long <- list(); ARG.short <- list(); ARG.cl = 1; ARG.cs = 1
MPH.long <- list(); MPH.short <- list(); MPH.cl = 1; MPH.cs = 1
AAR.long <- list(); AAR.short <- list(); AAR.cl = 1; AAR.cs = 1

seed.length <- c(ARG.long.dirs,   ARG.short.dirs,    ArgAla.long.dirs, ArgAla.short.dirs, 
		 		MetPhe.long.dirs, MetPhe.short.dirs, BES.long.dirs,    BES.short.dirs, 
		 		R5X.long.dirs,    R5X.short.dirs)

files <- c(ARG.long.dirs, ArgAla.long.dirs, MetPhe.long.dirs, BES.long.dirs, R5X.long.dirs)

REScentroid <- pdb_res_to_centroid("F:/3ebh.pdb")

REScentroid[,1] <- paste(REScentroid$resname, REScentroid$resnumber)
REScentroid[,1] <- gsub(" ", "\n", REScentroid[,1])
REScentroid     <- REScentroid[,c(1,3:5)]

ene         <- combine_energies(files)
energyRelations <- energy_relations(ene, files)
energyRelations$from <- as.vector(as.numeric(as.character(energyRelations$from)))
a <- 1:49
b <- 2:50
c <- rep(0,49)

new <- matrix(c(a,b,c),ncol=3)
new <- as.data.frame(new)
colnames(new)   <- colnames(energyRelations)
energyRelations <- rbind(energyRelations, new)

RESlength <- length(which(energyRelations$energy>=100))

quartileCOL 	<- quantile(energyRelations$energy[1:RESlength])
Low 	<- which(energyRelations$energy >=quartileCOL[1] & energyRelations$energy  <= quartileCOL[2])
Medium  <- which(energyRelations$energy > quartileCOL[2] & energyRelations$energy  <= quartileCOL[3])
High 	<- which(energyRelations$energy > quartileCOL[3] & energyRelations$energy  <= quartileCOL[4])
Highest <- which(energyRelations$energy > quartileCOL[4] & energyRelations$energy  <= quartileCOL[5])


lig.path  <- readCOM("F:/ENERGY_CALCULATIONS/Energy/HeatMap/new", pmat)
frame.len <- length(lig.path[[1]]$x)
interval.steps <- seq(1,frame.len,25)
inter.len <- length(interval.steps)-1

Xmean <- list()
Ymean <- list()
XmeanMin <- list()
YmeanMin <- list()
XmeanMax <- list()
YmeanMax <- list()

xtemp <- vector()
ytemp <- vector()
xtempMin <- vector()
ytempMin <- vector()
xtempMax <- vector()
ytempMax <- vector()

inc   = 25
from1 = 0
to1   = 0

for(s in 1:(length(lig.path)-1)){
	from1 = 0
	to1   = 0
	for(i in 1:inter.len){
		if(i == 1){
			from1 <- 1
			to1   <- inc + to1
			xtemp[i] <- mean(lig.path[[s]]$x[from1:to1])
			ytemp[i] <- mean(lig.path[[s]]$y[from1:to1])
			xtempMin[i] <- min(lig.path[[s]]$x[from1:to1])
			ytempMin[i] <- min(lig.path[[s]]$y[from1:to1])
			xtempMax[i] <- max(lig.path[[s]]$x[from1:to1])
			ytempMax[i] <- max(lig.path[[s]]$y[from1:to1])
		} else {
			from1  <- inc + from1
			to1	   <- inc + to1
			xtemp[i] <- mean(lig.path[[s]]$x[from1:to1])
			ytemp[i] <- mean(lig.path[[s]]$y[from1:to1])
			xtempMin[i] <- min(lig.path[[s]]$x[from1:to1])
			ytempMin[i] <- min(lig.path[[s]]$y[from1:to1])
			xtempMax[i] <- max(lig.path[[s]]$x[from1:to1])
			ytempMax[i] <- max(lig.path[[s]]$y[from1:to1])
		}	
	}
    Xmean[[s]] <- xtemp
    Ymean[[s]] <- ytemp
    XmeanMin[[s]] <- xtempMin
	YmeanMin[[s]] <- ytempMin
	XmeanMax[[s]] <- xtempMax
	YmeanMax[[s]] <- ytempMax
}

xmin.mat <- matrix(nrow=length(Xmean[[1]]) ,ncol=length(Xmean))
ymin.mat <- matrix(nrow=length(Xmean[[1]]) ,ncol=length(Xmean))
xmax.mat <- matrix(nrow=length(Xmean[[1]]) ,ncol=length(Xmean))
ymax.mat <- matrix(nrow=length(Xmean[[1]]) ,ncol=length(Xmean))

for(i in 1:length(Xmean[[1]])){
	for(j in 1:length(Xmean)){
		xmin.mat[i,j] <- XmeanMin[[j]][i]
		ymin.mat[i,j] <- YmeanMin[[j]][i]
		xmax.mat[i,j] <- XmeanMax[[j]][i]
		ymax.mat[i,j] <- YmeanMax[[j]][i]
	}
}
#
inter.len2 <- length(Xmean)
Xmean <- Reduce("+", Xmean); Xmean = Xmean/inter.len2
Ymean <- Reduce("+", Ymean); Ymean = Ymean/inter.len2
#
XmeanMin <- apply(xmin.mat,1,min)
YmeanMin <- apply(ymin.mat,1,min)
XmeanMax <- apply(xmax.mat,1,max)
YmeanMax <- apply(ymax.mat,1,max)

x2 <- as.matrix(REScentroid2[,2])
y2 <- as.matrix(REScentroid2[,3])
z2 <- as.matrix(REScentroid2[,4])

pamt 		<- matrix(0,nrow=4,ncol=4)
diag(pamt) 	<- 1:4
XY 			<- trans3D(x2, y2, z2, pmat = pamt) 



layoutLONG  <- as.matrix(XY[[1]]); layoutLONG <- cbind(layoutLONG, XY[[2]])
layoutLONG <- layoutLONG[indFORocc,]
colnames(layoutLONG) <- c("x", "y")
BesPullPath <- cbind(Xmean, Ymean); colnames(BesPullPath) <- colnames(layoutLONG)
layoutLONG  <- rbind(layoutLONG, BesPullPath)


#normalDist = (Distance-min(Distance))/(max(Distance)-min(Distance))

VSize <- vector()
ESize <- vector()
counter = 1
for(i in energyRelations$energy[1:RESlength])
{
	VSize[counter] <- (i-min(energyRelations$energy[1:RESlength]))/(max(energyRelations$energy[1:RESlength])-min(energyRelations$energy[1:RESlength]))
	counter = counter+1
}
VSize <- (VSize+1.8)^3.9; VSize <- VSize + 2
ESize <- as.vector(ListNodeOCC)
ESize <- (ListNodeOCC-min(ListNodeOCC))/(max(ListNodeOCC)-min(ListNodeOCC))
ESize <- (ESize+0.6)^3; ESize <- ESize ; ESize <- as.vector(ESize)
ESize <- ESize[indFORocc]
lohg <- length(which(energyRelations$energy>0))+1
LIGlength = 49
g <- graph.data.frame(energyRelations, directed=FALSE)
LbesPath <- length(V(g))-length(which(energyRelations$energy>=100))
V(g)$size=c(VSize, rep(2, LbesPath))
E(g)$size=c(ESize, rep(10, LbesPath-1))
VLabels <- c(REScentroid[indexing,1],rep("", LbesPath))
#
V(g)$color[Low]		<- "grey70"
V(g)$color[Medium]	<- "grey70"
V(g)$color[High]	<- "tan1"
V(g)$color[Highest]	<- "tomato3"
#
colfunc <- colorRampPalette(c("springgreen", "cyan", "blue", "darkblue"))
scaleF  <- colfunc(50)
#
colFUN <- V(g)$color
colFUN[lohg:(lohg+LbesPath)] <-  scaleF
colFUN <- colFUN[1:(length(colFUN)-1)]
E(g)$color = colFUN
V(g)$color = colFUN
V(g)$label.cex 	= 0.7
#
# POSITIVLY CHARGED RESIDUES #
t1 <- grep("LYS", VLabels)
t2 <- grep("ARG", VLabels)
pos.index <- c(t1, t2)
# NEGATIVLY CHARGED RESIDUES #
t1 <- grep("ASP", VLabels)
t2 <- grep("GLU", VLabels)
neg.index <- c(t1, t2)
# ACTIVE SITE RESIDUES #
ASres <- c("317", "319", "320","459" , "460", "461", "462", "463", "492", "493", "497", "518", "575", "580", "1034", "496", "500", "519", "526", "572")
as.index <- list()
for(i in 1:length(ASres))
{
	as.index[[i]] <- grep(ASres[i], VLabels)
}
as.index <- sort(unlist(as.index))

V(g)$framecolor[pos.index] <- "turquoise2"
V(g)$framecolor[neg.index] <- "darkred"
V(g)$framecolor[as.index]  <- "black"
V(g)$framecolor[(lohg+1):(lohg:lohg+50)] <- NA
i.i.i <- which(is.na(V(g)$framecolor))
V(g)$framecolor[i.i.i] <- ""
# Plot Network
id2 <- tkplot(g, layout=layoutLONG, rescale=TRUE, 
	vertex.label = VLabels, 
	vertex.shape = "fcircle",
	vertex.label.cex=0.5,
	vertex.label.color="black",
	edge.width=E(g)$size,
	vertex.frame.color=V(g)$framecolor,
	edge.color=E(g)$color,
	edge.curved=seq(-0.2, 0.2))











