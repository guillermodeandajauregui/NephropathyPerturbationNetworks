genelevel_xtalk_network<-function(list_graphs, pathways){
  #a list of igraph gene level pathway graphs, NAMED
  #a list of pathway names of interest (from a crosstalk network)
  #names should match (else, use fuzzy matching to generate fuzzy names)
  
  #subset list of pathway graphs
  list_graphs = list_graphs[names(list_graphs)%in%pathways]
  
  #merge the first two pathways
  xtalk_nw = igraph::union(list_graphs[[1]], list_graphs[[2]])
  
  #merge the rest
  for(i in 3:length(list_graphs)){
    xtalk_nw = igraph::union(xtalk_nw, list_graphs[[i]])
  }
  rm(i)
  
  #simplify network
  xtalk_nw = igraph::simplify(graph = xtalk_nw, remove.multiple = TRUE, remove.loops = TRUE)
  return(xtalk_nw)
}
