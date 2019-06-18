#Installing the package
install.packages("igraph")

#Data: Protein-Protein Interaction Network
library(igraph)
biogrid <- read.delim("data/BIOGRID.txt", stringsAsFactors = F)    

#Look the structure
str(biogrid)
names(biogrid)
gene.data <- data.frame(biogrid$Entrez.Gene.Interactor.A, biogrid$Entrez.Gene.Interactor.B)
HSnet <- graph.data.frame(gene.data, directed=F)

#plotting the network
plot(HSnet)

#Adjacency of the graph network
A <- get.adjacency(HSnet)
str(A)
length(A)

#Multiple edges
A[1:15,1:15]

# the following is FALSE if the graph is not simple
is.simple(HSnet)

# remove multiple edges and self-loops
HSnet <- simplify(HSnet, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = getIgraphOpt("edge.attr.comb"))
is.simple(HSnet)
A <- get.adjacency(HSnet)

# only single edges now
A[1:15,1:15]
length(A)

# for this application, we remove nodes of very high degree; these are usually house-keeping
#genes that are necessary to keep a cell alive, but are usually not specific to a particular disease.
overly.attached.proteins <- which(degree(HSnet)>1000) 
HSnet <- delete.vertices(HSnet, overly.attached.proteins)

# the following is TRUE if the graph is connected.
is.connected(HSnet)

