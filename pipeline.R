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

##pathway annotation
load(inputs[2])


a = ReGAGE(expmatrix = matriz, 
           pathways = PWs[1:20], 
           cases = dbdb_Glom_num, 
           controls = Cont_Glom_num, 
           qvalue = 0.1, 
           setsize = 10)
