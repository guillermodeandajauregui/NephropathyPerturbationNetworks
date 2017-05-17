##########################################
#Pipeline for the exploration of 
#Pathway deregulation in the 
#Diabetic Nephropathy model
##########################################

#inputs 

## expression matrix 
## grouping information
## pathways

#outputs 

## Pathway Perturbation Networks + communities 
## Driving genes 

#1) Pathway and crosstalk perturbation > ReGAGEX
#2) Pw Pert NW analysis > NetworkAnalyzer Clone 
#3) Null model generation
#4) Community evolution  > community evolution functions
#5) Gene-level crosstalk network > Fuzzy Matcher + Graphite + Igraph
#6) Driver identification > DeMAND

#sources and libraries
inputs = scan(file = "inputs.txt", what = "character")
source(file = "sources.R")

# load inputs 

##use expression matrix with older normalization, for consistency with the neuropathy results 
matriz = fread(input = inputs[1], data.table = FALSE)
rownames(matriz) <- matriz[,1]
matriz <- matriz[,-1]
matriz <- matriz[order(rownames(matriz)),]
load(inputs[4])
matriz <- matriz[stringr::str_to_upper(rownames(matriz))%in%features_exp_norm,] #remove non protein coding genes

##groups
groups = fread(inputs[3])
table(groups$Group)

Cont_Glom_num = which(colnames(matriz)%in%groups[Group == "Cont_Glom"]$Sample)
dbdb_Glom_num = which(colnames(matriz)%in%groups[Group == "dbdb_Glom"]$Sample)
dbdb_Pio_Glom_num = which(colnames(matriz)%in%groups[Group == "dbdb-Pio_Glom"]$Sample)

Cont_Cortex_num = which(colnames(matriz)%in%groups[Group == "Cont_Cortex"]$Sample)
dbdb_Cortex_num = which(colnames(matriz)%in%groups[Group == "dbdb_Cortex"]$Sample)
dbdb_Pio_Cortex_num = which(colnames(matriz)%in%groups[Group == "dbdb-Pio_Cortex"]$Sample)

Cont_SCN_num = which(colnames(matriz)%in%groups[Group == "Cont_SCN"]$Sample)
dbdb_SCN_num = which(colnames(matriz)%in%groups[Group == "dbdb_SCN"]$Sample)
dbdb_Pio_SCN_num = which(colnames(matriz)%in%groups[Group == "dbdb-Pio_SCN"]$Sample)

##pathway annotation
#load(inputs[2])
reactome_pw = pathways(species = "mmusculus", "reactome")
reactome_pw = convertIdentifiers(x = reactome_pw, to = "symbol")
reactome_pw = lapply(reactome_pw, function(x) nodes(pathwayGraph(x)))


#copy matriz 
#matriz_copia = matriz
#rownames(matriz_copia) <- stringr::str_to_upper(rownames(matriz_copia))
#a = gage(exprs = matriz, gsets = PWs[1:20], ref = Cont_Glom_num, samp = dbdb_Glom_num, set.size = 10, same.dir = TRUE, compare = "unpaired")
#a = ReGAGE(expmatrix = matriz_copia, pathways = PWs, cases = dbdb_Glom_num, controls = Cont_Glom_num, qvalue = 0.1)
# components(a$g)
# b = ReGAGE(expmatrix = matriz, pathways = reactome_pw, cases = dbdb_SCN_num, controls = Cont_SCN_num, qvalue = 0.05)
# plot(b$g)
# components(b$g)
# V(b$g)$name

###############################################################################
#1) Pathway and crosstalk perturbation > ReGAGEX
###############################################################################

#ReGAGE, qvalue 0.05


#glom
glom_regages = list(alfa = ReGAGE(expmatrix = matriz, 
                                  pathways = reactome_pw, 
                                  cases = dbdb_Glom_num, 
                                  controls = Cont_Glom_num, 
                                  qvalue = 0.05),
                    beta = ReGAGE(expmatrix = matriz, 
                                  pathways = reactome_pw, 
                                  cases = dbdb_Pio_Glom_num, 
                                  controls = dbdb_Glom_num, 
                                  qvalue = 0.05),
                    delta = ReGAGE(expmatrix = matriz, 
                                  pathways = reactome_pw, 
                                  cases = dbdb_Pio_Glom_num, 
                                  controls = Cont_Glom_num, 
                                  qvalue = 0.05)
                    )

#cortex
Cortex_regages = list(alfa = ReGAGE(expmatrix = matriz, 
                                  pathways = reactome_pw, 
                                  cases = dbdb_Cortex_num, 
                                  controls = Cont_Cortex_num, 
                                  qvalue = 0.05),
                    beta = ReGAGE(expmatrix = matriz, 
                                  pathways = reactome_pw, 
                                  cases = dbdb_Pio_Cortex_num, 
                                  controls = dbdb_Cortex_num, 
                                  qvalue = 0.05),
                    delta = ReGAGE(expmatrix = matriz, 
                                   pathways = reactome_pw, 
                                   cases = dbdb_Pio_Cortex_num, 
                                   controls = Cont_Cortex_num, 
                                   qvalue = 0.05)
)

#SCN - new version
SCN_regages = list(alfa = ReGAGE(expmatrix = matriz, 
                                    pathways = reactome_pw, 
                                    cases = dbdb_SCN_num, 
                                    controls = Cont_SCN_num, 
                                    qvalue = 0.05),
                      beta = ReGAGE(expmatrix = matriz, 
                                    pathways = reactome_pw, 
                                    cases = dbdb_Pio_SCN_num, 
                                    controls = dbdb_SCN_num, 
                                    qvalue = 0.05),
                      delta = ReGAGE(expmatrix = matriz, 
                                     pathways = reactome_pw, 
                                     cases = dbdb_Pio_SCN_num, 
                                     controls = Cont_SCN_num, 
                                     qvalue = 0.05)
)

#save(SCN_regages, Cortex_regages, glom_regages, file = "...PIO/regages.RData")
#save.image("...PIO/workspace_phase1.RData")


###############################################################################
#2) Pw Pert NW analysis > NetworkAnalyzer Clone 
###############################################################################

glom_analysis = lapply(glom_regages, FUN = function(x){
   NetworkAnalyzer(x$g, 
                   directed = FALSE, 
                   skip.betweenness = FALSE, 
                   workaround.betweenness = TRUE)
 })

Cortex_analysis = lapply(Cortex_regages, FUN = function(x){
  NetworkAnalyzer(x$g, 
                  directed = FALSE, 
                  skip.betweenness = FALSE, 
                  workaround.betweenness = TRUE)
})

SCN_analysis = lapply(SCN_regages, FUN = function(x){
  NetworkAnalyzer(x$g, 
                  directed = FALSE, 
                  skip.betweenness = FALSE, 
                  workaround.betweenness = TRUE)
})

#save(SCN_analysis, Cortex_analysis, glom_analysis, file = "...PIO/nw_analysis_regage.RData")

glom_nws = lapply(glom_analysis, function(x){
  x$g
})

Cortex_nws = lapply(Cortex_analysis, function(x){
  x$g
})

SCN_nws = lapply(SCN_analysis, function(x){
  x$g
})

#save(SCN_nws, Cortex_nws, glom_nws, file = "...PIO/regage_analyzed_nw_list.RData")
#save.image("...PIO/workspace_phase2.RData")

###############################################################################
#3) Null model generation
###############################################################################

#background probability crosstalk network

background_network_reactome = jaccard_matrix(reactome_pw, reactome_pw)
background_network_reactome = graph_from_adjacency_matrix(background_network_reactome, weighted = TRUE, mode = "undirected")
background_network_reactome = simplify(background_network_reactome, remove.multiple = TRUE, remove.loops = TRUE)
background_network_reactome = NetworkAnalyzer(background_network_reactome)
distribution_plot(background_network_reactome$g)

glom_bootstraps = lapply(X = glom_nws, FUN = function(g){
  bootstrap_model_4_function(g = g, background = background_network_reactome$g, n = 5000)
})

cortex_bootstraps = lapply(X = Cortex_nws, FUN = function(g){
  bootstrap_model_4_function(g = g, background = background_network_reactome$g, n = 5000)
})

SCN_bootstraps = lapply(X = SCN_nws, FUN = function(g){
  bootstrap_model_4_function(g = g, background = background_network_reactome$g, n = 5000)
})

#save(glom_bootstraps, cortex_bootstraps, SCN_bootstraps, background_network_reactome, file = "...PIO/background_and_bootstraps.R")
#save.image("...PIO/workspace_phase3.RData")

#### bootstraps with a non-normalized version of the Jaccard probability; they don't diverge wildly so unnecessary
# glom_bootstraps_nonnormalized = lapply(X = glom_nws, FUN = function(g){
#   bootstrap_model_4_function(g = g, background = background_network_reactome$g, n = 5000)
# })
# 
# cortex_bootstraps_nonnormalized = lapply(X = Cortex_nws, FUN = function(g){
#   bootstrap_model_4_function(g = g, background = background_network_reactome$g, n = 5000)
# })
# 
# SCN_bootstraps_nonnormalized = lapply(X = SCN_nws, FUN = function(g){
#   bootstrap_model_4_function(g = g, background = background_network_reactome$g, n = 5000)
# })


# comparative_degree_distribution_plot_points(nw = glom_nws$delta, bootstrap = glom_bootstraps$delta)
# comparative_degree_distribution_plot_points(nw = Cortex_nws$delta, bootstrap = cortex_bootstraps$delta)
# comparative_degree_distribution_plot_points(nw = SCN_nws$delta, bootstrap = SCN_bootstraps$delta)

# for(i in 1:length(glom_bootstraps$alfa)){
# plot_nicely(NetworkAnalyzer(glom_bootstraps$alfa[[i]], directed = FALSE, skip.betweenness = FALSE, workaround.betweenness = TRUE))
# }
# plot_nicely(glom_analysis$alfa)
# 
# smallworld_test(g = glom_analysis$alfa$g, NullModelList = glom_bootstraps$alfa)

###############################################################################
#4) Community evolution  > community evolution functions
###############################################################################

#alfa to delta ... how connectivity changes from the disease state to the treated state
#beta to delta ... differences in connectivity of alterations inside the disease state, against disease v non disease

#alfa to delta ... how connectivity changes from the disease state to the treated state
glom_EvoComm_alfa_delta = network_community_evolution_analysis(gi = glom_nws$alfa, 
                                                               gj = glom_nws$delta, 
                                                               k = 0.2, 
                                                               grouping = "infomap")

Cortex_EvoComm_alfa_delta = network_community_evolution_analysis(gi = Cortex_nws$alfa, 
                                                               gj = Cortex_nws$delta, 
                                                               k = 0.2, 
                                                               grouping = "infomap")

SCN_EvoComm_alfa_delta = network_community_evolution_analysis(gi = SCN_nws$alfa, 
                                                               gj = SCN_nws$delta, 
                                                               k = 0.2, 
                                                               grouping = "infomap")


#beta to delta ... differences in connectivity of alterations inside the disease state, against disease v non disease
glom_EvoComm_beta_delta = network_community_evolution_analysis(gi = glom_nws$beta, 
                                                               gj = glom_nws$delta, 
                                                               k = 0.2, 
                                                               grouping = "infomap")

Cortex_EvoComm_beta_delta = network_community_evolution_analysis(gi = Cortex_nws$beta, 
                                                                 gj = Cortex_nws$delta, 
                                                                 k = 0.2, 
                                                                 grouping = "infomap")

SCN_EvoComm_beta_delta = network_community_evolution_analysis(gi = SCN_nws$beta, 
                                                              gj = SCN_nws$delta, 
                                                              k = 0.2, 
                                                              grouping = "infomap")

SCN_EvoComm_beta_delta$comms_gj

# save(glom_EvoComm_alfa_delta, 
#      glom_EvoComm_beta_delta, 
#      Cortex_EvoComm_alfa_delta, 
#      Cortex_EvoComm_beta_delta, 
#      SCN_EvoComm_alfa_delta, 
#      SCN_EvoComm_beta_delta, 
#      file = "...PIO/EvoComms.RData")


###############################################################################
#5) Gene-level crosstalk network > Fuzzy Matcher + Graphite + Igraph
###############################################################################

reactome_graphs = pathways(species = "mmusculus", "reactome")
reactome_graphs = convertIdentifiers(x = reactome_graphs, to = "symbol")
reactome_graphs = lapply(reactome_graphs, pathwayGraph)
reactome_graphs = lapply(reactome_graphs, igraph::graph_from_graphnel)

glom_gene_xtalk_nw = lapply(X = glom_nws, function(g){
  genelevel_xtalk_network(list_graphs = reactome_graphs, pathways = names(V(g)))
})

Cortex_gene_xtalk_nw = lapply(X = Cortex_nws, function(g){
  genelevel_xtalk_network(list_graphs = reactome_graphs, pathways = names(V(g)))
})

SCN_gene_xtalk_nw = lapply(X = SCN_nws, function(g){
  genelevel_xtalk_network(list_graphs = reactome_graphs, pathways = names(V(g)))
})

###############################################################################
#6) Driver identification > DeMAND
###############################################################################

#for each tissue, for each comparison (9 networks total)

####################################
############make "annotation" file
####################################


demand_anot = cbind(rownames(matriz), rownames(matriz))

####################################
############make demand lists
####################################
glom_demand =         lapply(glom_gene_xtalk_nw, function(x){
                      demandClass(exp = matriz, 
                                  anno = demand_anot, 
                                  network = get.edgelist(x)
                      )
                    }
                    )

Cortex_demand =         lapply(Cortex_gene_xtalk_nw, function(x){
  demandClass(exp = matriz, 
              anno = demand_anot, 
              network = get.edgelist(x)
  )
}
)

SCN_demand =         lapply(SCN_gene_xtalk_nw, function(x){
  demandClass(exp = matriz, 
              anno = demand_anot, 
              network = get.edgelist(x)
  )
}
)

####################################
############make lists of fg / bg
####################################


glom_contrasts = list(
  alfa = list(fg = dbdb_Glom_num, bg = Cont_Glom_num),
  beta = list(fg = dbdb_Pio_Glom_num, bg = dbdb_Glom_num),
  gamma = list(fg = dbdb_Glom_num, bg = Cont_Glom_num)
)

Cortex_contrasts = list(
  alfa = list(fg = dbdb_Cortex_num, bg = Cont_Cortex_num),
  beta = list(fg = dbdb_Pio_Cortex_num, bg = dbdb_Cortex_num),
  gamma = list(fg = dbdb_Cortex_num, bg = Cont_Cortex_num)
)

SCN_contrasts = list(
  alfa = list(fg = dbdb_SCN_num, bg = Cont_SCN_num),
  beta = list(fg = dbdb_Pio_SCN_num, bg = dbdb_SCN_num),
  gamma = list(fg = dbdb_SCN_num, bg = Cont_SCN_num)
)

# glom_demand_run = mapply(FUN = function(demanda, contrastes){
#   runDEMAND2(x = demanda, 
#              fgIndex = contrastes[[1]], 
#              bgIndex = contrastes[[2]], 
#              method = "integers")
#             },
#   glom_demand,
#   glom_contrasts
#                 )

#########################################
############Run 10 DeMANDS per comparison
#########################################


glom_demand_bootstrap = mapply(FUN = function(demanda, contrastes){
  x = 1:2
  lapply(X = x, FUN = function(x){
    runDEMAND2(x = demanda, 
               fgIndex = contrastes[[1]], 
               bgIndex = contrastes[[2]], 
               method = "integers")
  })
},
demanda = glom_demand,
contrastes = glom_contrasts, 
SIMPLIFY = FALSE
)

Cortex_demand_bootstrap = mapply(FUN = function(demanda, contrastes){
  x = 1:2
  lapply(X = x, FUN = function(x){
    runDEMAND2(x = demanda, 
               fgIndex = contrastes[[1]], 
               bgIndex = contrastes[[2]], 
               method = "integers")
  })
},
demanda = Cortex_demand,
contrastes = Cortex_contrasts, 
SIMPLIFY = FALSE
)

SCN_demand_bootstrap = mapply(FUN = function(demanda, contrastes){
  x = 1:2
  lapply(X = x, FUN = function(x){
    runDEMAND2(x = demanda, 
               fgIndex = contrastes[[1]], 
               bgIndex = contrastes[[2]], 
               method = "integers")
  })
},
demanda = SCN_demand,
contrastes = SCN_contrasts, 
SIMPLIFY = FALSE
)

save.image("~/X/PIO/workspace_phase6.RData")

#########################################
############To explore results
#########################################

Reduce(intersect, lapply(SCN_demand_bootstrap$beta, function(x){
  x@moa$moaGene[1:50]
}))

#########################################
############DONE!########################
#########################################