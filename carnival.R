# In this code we will go through the network propagation in the UNITI or in the UNIFI dataset
# The end goal would be to develop a shiny R app from it.
# Input normalised transcriptomic data either DESEQ2 vsn or RMA normed microarrays
# Necessary packages:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("OmnipathR")
BiocManager::install("dorothea")
BiocManager::install("CARNIVAL")
library("OmnipathR")
library("dorothea")
library("viper")
library(decoupleR)
library(igraph)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(CARNIVAL)

file.exists("D:/RProject/Project/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe")

# Reading in Flod Change data
DE_data <- read.csv("D:/PythonProject/CIProject/resource/TCGA/DEA_FC_pvalue.csv")
DE_data$FDR <- NULL
### ERROR: some DEGs have multiple FC values (probably due to isoforms) 
### -- SOLUTION: calculate the average
DE_data <- DE_data %>%
  group_by(gene_name) %>%
  summarise(average_value = mean(logFC))
DE_data <- data.frame(row.names = DE_data$gene_name, logFC = DE_data$average_value)
colnames(DE_data)[colnames(DE_data)== "logFC"] <- "FC"

# Calculating transcription factor weights specific for that particular patient 
# Downloading regulatory network
net_collectri <- decoupleR::get_collectri(split_complexes = FALSE)
net_collectri
?decoupleR::get_collectri
#Running a generalised linear model on the transcriptomic data using 
# the decoupler package 
# The run_ulm is a linear model, it replaces the VIPER analysis
#Calculate TF activity
?run_ulm
# TF has at least 5 interactions
contrast_acts  <- run_ulm(mat=DE_data, net=net_collectri, 
                       .source='source', .target='target',
                       .mor='mor', minsize = 5)
?run_ulm
# keep significantly different activity between the condition (Normal vs. Tumour)
contrast_acts <- contrast_acts[contrast_acts$p_value <0.05, ]
contrast_acts

#We will use the tibble f_contrast_acts as the Carnival input.

# 6. Filtering the protein interaction network

# The interaction network needs to be filtered by expression values.
# Creating an OmniPath PPI interaction network. We are using all the nodes as 
# they are we do not filter out complexes here first. We filter directed and 
# singed interactions and define interaction weight as stimulation - inhibition.
# Meaning that inhibition has a unit negative weight. 
interaction_network <- import_omnipath_interactions(organism = 9606)
interaction_network
interaction_network <- interaction_network[(interaction_network$consensus_direction == 1) & 
                                             ((interaction_network$consensus_stimulation == 1) | 
                                                (interaction_network$consensus_inhibition == 1)), 
                                           c("source_genesymbol", "target_genesymbol","consensus_stimulation", "consensus_inhibition","n_references")]

interaction_network$wheight <- interaction_network$consensus_stimulation - interaction_network$consensus_inhibition 

interaction_network_ready <- data.frame("source" = interaction_network$source_genesymbol,
                                        "target"= interaction_network$target_genesymbol,
                                          "interaction" = interaction_network$wheight)
interaction_network_ready <- interaction_network_ready[interaction_network_ready$interaction != 0,]

# READING THE TUMOUR CONDITION EXPRESSED GENES
expression <- read.csv("D:/PythonProject/CIProject/resource/TCGA/LogAveTP_genesymbol.csv")
#The selection cut off will be the mean - 2 SD. This means that we have all the
# expressed genes which have even a faint signal in our microarray dataset. 
#Also the selected state is the control state. That is what we target with drug.
# calculate the z score!
expression
mean_expr =  mean(expression$Log2.x.1.)
mean_expr
sd_expr <- sd(expression$Log2.x.1, na.rm = TRUE)
sd_expr
expressed_gene <- expression[expression$Log2.x.1. > (mean_expr - 2 * sd_expr),]
expressed_gene <- expressed_gene[expressed_gene$Log2.x.1.!=0,]
expressed_gene
#expressed_gene <- expressed_gene[order(expressed_gene$Log2.x.1., decreasing = TRUE),]
#percent <- ceiling(nrow(expressed_gene)*0.75)
#expressed_gene_75 <- expressed_gene[1:percent,]

#Here we filter to complexes, because gene expression is gene based and not complex based we do not work with complexes. 
interaction_network_specific <- interaction_network_ready[interaction_network_ready$source %in% expressed_gene$Gene, ] 
interaction_network_specific <- interaction_network_specific[interaction_network_specific$target %in% expressed_gene$Gene, ]
interaction_network_ready
express_and_in_interaction <- interaction_network_ready[interaction_network_ready$source %in% expressed_gene$Gene, "source"]
express_and_in_interaction <- union(express_and_in_interaction,interaction_network_ready[interaction_network_ready$target %in% expressed_gene$Gene, "target"])
express_and_in_interaction

contrast_acts_f <- contrast_acts[contrast_acts$source %in% express_and_in_interaction,]
contrast_acts_f
transcriptional_signal <- contrast_acts_f$score

names(transcriptional_signal) <- contrast_acts_f$source

#Using only the giant component of the expressed graph.
g<- graph_from_data_frame(interaction_network_specific)
components <- igraph::components(g, mode="weak")
biggest_cluster_id <- which.max(components$csize)

# ids
vert_ids <- V(g)[components$membership == biggest_cluster_id]

# subgraph
g_giant_component <- igraph::induced_subgraph(g, vert_ids)
interaction_network_giant <- igraph::as_data_frame(g_giant_component, what = "edges")
interaction_network_giant$source <- interaction_network_giant$from
interaction_network_giant$interaction <- interaction_network_giant$interaction
interaction_network_giant$target <- interaction_network_giant$to
interaction_network_giant <-  subset(interaction_network_giant, 
                                     select = c("source", "interaction","target"))

# 7. Creating targets and running CARNIVAL on the network
install.packages("mygene")
if(!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("mygene")
library(mygene)

#perturbation point
iupred <- read.csv("D:/PythonProject/CIProject/result/iupred_output.csv", sep = '\t')
uniprot_ids <- iupred$X..Human.Protein

target <- queryMany(uniprot_ids, scopes='uniprot', fields="symbol", species="human") 
target <- unique(target$symbol)
vec = rep(1, length(target))
names(vec) <- target
vec


# Print the results
# print(target) 
?defaultCplexCarnivalOptions
carnivalOptions <- defaultCplexCarnivalOptions()
carnivalOptions$solverPath <- "D:/RProject/Project/CPLEX_Studio2211/cplex/bin/x64_win64/cplex.exe"
carnivalOptions$threads <- 20
carnivalOptions$timelimit <- 3600*0.5
carnivalOptions$keepLPFiles <- FALSE #keep the file or not
carnivalOptions$solver <- "cplex" 
carnivalOptions$cplexMemoryLimit <- 32000
carnivalOptions$workdir <- "D:/RProject/Project/CARNIVAL"
carnivalOptions$outputFolder <- "D:/RProject/Project/CARNIVAL"
?carnivalOptions$solver
checkOptionsValidity(solver = "cplex")

#RUN
resultsLpSolve <- runVanillaCarnival( perturbations = vec,
                                      measurements = transcriptional_signal,
                                      priorKnowledgeNetwork = interaction_network_giant, 
                                      carnivalOptions = carnivalOptions)
?runVanillaCarnival
write.table(resultsLpSolve$nodesAttributes, "D:/RProject/Project/CARNIVAL/Results/CARNIVAL_node_attribute.txt", sep= "\t" , quote= FALSE)
write.table(resultsLpSolve$weightedSIF, "D:/RProject/Project/CARNIVAL/Results/CARNIVALnetwork.txt", sep= "\t" , quote= FALSE)


#To do the proper analysis we need to choose the samples:
#metadata_selected <- metadata[(metadata$TRT01A != "Placebo IV") &  (metadata$TRT01A != "Normal") & (metadata$location=="ILEUM"),]
