#Data: Genes causing autism
genes.table <- read.delim("data/gene-id-table.txt")   
names(genes.table)

# read the scores for Autism
gene.score<-read.csv("data/gene-score.csv",stringsAsFactors=F) 
attach(gene.score)
names(gene.score)

#Display unique scores
unique(Score)

#Identifying genes that have significant scores
signif.scores<-c("3","1S","1","2S","2","3S")
signif.genes<-Gene.Symbol[which(Score %in% signif.scores)] 
signif.EIDs <- genes.table[which(genes.table[,1] %in% signif.genes),2]

# Now use the protein interaction network HSnet, created previously, 
#to determine the genes that are present in the network and known to cause Autism
geneEIDs <- as.numeric(V(HSnet)$name)
HSnetN<-HSnet
V(HSnetN)$name<-1:length(V(HSnet))
signif.ids<-which(geneEIDs %in% signif.EIDs) 
length(signif.ids)

#Building an Autism Interactome
source('steiner_tree.R')

# Identify the Steiner Tree and note the time this function call takes
system.time(HS.stree <- steiner_tree(terminals=signif.ids, graph=HSnetN))

# identify the genes that have significant scores, and assign the color “red” to them
colors<-rep("skyblue",length(V(HS.stree)))
colors[which(as.numeric(V(HS.stree)$name) %in% signif.ids)] = "red"

# assign colors to the vertices of the tree
V(HS.stree)$color = colors

# plot and save to file
pdf("ASD_interactome.pdf",width=12, height=12) 
system.time(plot(HS.stree,vertex.label=labels,vertex.size=5,vertex.label.cex=0.8)) 
dev.off()