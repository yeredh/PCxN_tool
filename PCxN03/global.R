library(data.table)
path.net.dframe <- data.table(readRDS("data/pcxn_gobp.RDS"))
path.net.dframe$PathCor <- round(path.net.dframe$PathCor,4)
path.net.dframe$Overlap.Coefficient <- round(path.net.dframe$Overlap.Coefficient,4)


# Get the names for all pathways
path.names = unique(c(path.net.dframe$Pathway.A,path.net.dframe$Pathway.B))
path.names <- sort(path.names)
path.list <- as.list(path.names)
names(path.list) <- path.names