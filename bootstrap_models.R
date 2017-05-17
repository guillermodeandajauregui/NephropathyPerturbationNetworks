###########################################
#ReGAGE null model generation 
###########################################

#background network
## A network of pathways connected if they share molecules
## Edges are weighted (by Jaccard index between two pathways)

# jaccard_matrix(list_pathways, list_pathways)


#Bootstrap model
#takes a pathway perturbation network
#a background network, with all the pathways tested for the original enrichment analysis
#a number of repetitions

bootstrap_model_4_function <-function(g, background, n){
  bs = lapply(1:n, function(x){
    gi = induced_subgraph(graph = background, #keeps only the nodes in the perturbed network
                          vids = V(g)$name)
    gi = delete.edges(graph = gi, edges = sample(x = E(gi), 
                                                 size = length(E(gi)) - length(E(g)),
                                                 replace = FALSE, 
                                                 prob = 1 - (E(gi)$weight/max(E(gi)$weight))
                                                 #prob = 1 - E(gi)$weight #non normalized .. it doesn't change things brutally
                                                ) 
    )
    
    return(gi)
  }
  )
}

#use bootstrap model to evaluate whether graph has small-world properties
  #heuristic L/L_null ~ 1 C/C_null > 1
smallworld_test <- function(g, NullModelList){
  L_null = mean(sapply(NullModelList, mean_distance, directed = FALSE, unconnected = TRUE))
  L = mean_distance(g)
  
  C_null = mean(sapply(NullModelList, transitivity, type = "undirected"))
  C = transitivity(g, type = "undirected")
  
  return(list(l = L/L_null, c = C/C_null))
  
}

#function to plot degree distribution

distribution_plot = function(nw, coloir = "tomato2"){
   plot(x = names(table(degree(graph = nw))),
        y = table(degree(graph = nw)),
        col = coloir
  )
}

distro_plot = function(nw, coloir = "tomato2"){
  lines(x = names(table(degree(graph = nw))),
        y = table(degree(graph = nw)),
        col = coloir
  )
}

distro_plot_points = function(nw, coloir = "tomato2"){
  points(x = names(table(degree(graph = nw))),
        y = table(degree(graph = nw)),
        col = coloir
  )
}

comparative_degree_distribution_plot <- function(nw, bootstrap, n = length(bootstrap)){
  #makes a plot with the degree distribution of a network
  #along with the degree distributions of the bootstrap model
  
  #first, plot for network
  plot(x = names(table(degree(graph = nw))),
       y = table(degree(graph = nw)),
       xlim = c(0,max(as.numeric(names(table(degree(graph = nw)))))),
       #xaxp = c(0, 20, 20),
       ylim = c(0, max(table(degree(graph = nw)))),
       col = "tomato2",
       xlab = "degree",
       ylab = "frequency"
      )
  
  #make a palette of colors 
  colorz = sample(rainbow(n = n), n)
  k = 0
  for(i in bootstrap){
    k = k +1
    distro_plot(i, colorz[k])
  }
  
}

comparative_degree_distribution_plot_points <- function(nw, bootstrap, n = length(bootstrap)){
  #makes a plot with the degree distribution of a network
  #along with the degree distributions of the bootstrap model
  
  #first, plot for network
  plot(x = names(table(degree(graph = nw))),
       y = table(degree(graph = nw)),  type = "b",
       xlim = c(0,max(as.numeric(names(table(degree(graph = nw)))))),
       #xaxp = c(0, 20, 20),
       ylim = c(0, max(table(degree(graph = nw)))),
       col = "tomato2",
       xlab = "degree",
       ylab = "frequency"
  )
  
  #make a palette of colors 
  colorz = sample(rainbow(n = n), n)
  k = 0
  for(i in bootstrap){
    k = k +1
    distro_plot_points(i, colorz[k])
  }
  
}


