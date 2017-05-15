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

