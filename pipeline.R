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
#3) Gene-level crosstalk network > Fuzzy Matcher + Graphite + Igraph
#4) Driver identification > DeMAND

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
