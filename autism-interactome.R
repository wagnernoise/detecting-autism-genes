install.packages('sna')
library(sna, quietly=TRUE)


#A function that computes the connectivity scores for a network
#Aere the scores are diameters of the network and average geodesic distance between any two nodes
c.scores<-function(graph) {
  n<-length(V(graph))
  sp<-shortest.paths(graph)
  neighbors<-sum(sp==1)/2
  neighbors2<-sum(sp==2)/2 
  return(c(2*neighbors/(n*(n-1)),2*neighbors2/(n*(n-1))))
}

clus<-clusters(HSnetN, mode=c("weak")) 
connected.ids<-which(clus$membership==1) 
length(connected.ids)

#Generate N randomly chosen subnetworks. Note: this will take a while if N is set large.
N<-50 
strees<-list(N) 
effs<-numeric(N) 
nei<-numeric(N) 
nei2<-numeric(N) 

for (i in 1:N) {
  new.ids<-sample(x=connected.ids,size=length(signif.ids)) 
  strees[[i]] <- steiner_tree(terminals=new.ids, graph=HSnetN) 
  effs[i]<-efficiency(get.adjacency(strees[[i]],sparse=F)) 
  cs<-c.scores(strees[[i]])
  nei[i]<-cs[1]
  nei2[i]<-cs[2] 
}
# print the efficiencies and connectivity scores for each of the N random graphs
effs 
nei
nei2

# Finally, print the efficiency score and connectivity scores for the Autism Interactome
efficiency(get.adjacency(HS.stree,sparse=F)) 
c.scores(HS.stree)

#Identifying New Candidate Genes
# compute the betweeness centrality scores for each node
betweeness_centrality_scores = igraph::betweenness(HS.stree)

# now identify only those NOT already known to be significant
significant_centrality = c()
count = 0
for (i in 1:length(betweeness_centrality_scores)) {
  if (!(as.numeric(names(betweeness_centrality_scores[i])) %in% signif.ids)) { 
    significant_centrality = c(significant_centrality, betweeness_centrality_scores[i]) 
    }
}
# sort
significant_centrality = sort(significant_centrality, decreasing=TRUE) 
significant_centrality
length(significant_centrality)

