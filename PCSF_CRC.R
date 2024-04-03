library(OmnipathR)
library(PCSF)
library(topGO)
?topGO
##WEIGHT OMNIPATH NETWORK BASED ON NUMBER OF RESOURCES PROVING THE INTERACTION##
op_ppi <- import_omnipath_interactions()
resources <- op_ppi$n_resources

# Inverse transformation, interaction more reliable -> less cost
inverse_values <- 1 / resources

# Min-max normalization
min_val <- min(inverse_values)
max_val <- max(inverse_values)
normalized_values <- (inverse_values - min_val) / (max_val - min_val)

# Assign normalized values as cost to the edges
cost <- normalized_values

op_ppi$cost <- normalized_values

# Specify the columns you want to keep
columns_to_keep <- c("source_genesymbol", "target_genesymbol", "cost")

# Subset the dataframe to keep only the specified columns
op_ppi <- op_ppi[, columns_to_keep]


# Specify the new column names
new_column_names <- c("from", "to","cost")

# Rename the columns
colnames(op_ppi) <- new_column_names

op_ppi <- as.data.frame(op_ppi)

#construct an interaction network
ppi <- construct_interactome(op_ppi)
?construct_interactome


###########CRC data analysis####################################################

data <- read.csv('D:/PythonProject/CIProject/resource/TCGA/COAD_DEGs_noCF.csv')

#$Expression = name of column where the LogFC are, $GeneSymbol = name of column where the gene symbols are
terminals_crc <- setNames(data$logFC, data$gene_name)

#PCSF Analysis - return a igraph 'subnet'
subnet <- PCSF(ppi, terminals_crc, w = 2, b = 1, mu = 0.0005)
?PCSF
#Plot the subnetwork
plot.PCSF(subnet)
?igraph 
#Enrichment Analysis
res_enrichR <- enrichment_analysis(subnet)
?enrichment_analysis
##Plot an interactive subnetwork with functional enrichment analysis
plot.PCSFe(res_enrichR$subnet)
?plot.PCSFe

##############################WRITING OUT RESULTS##############################

node_attributes <- data.frame(name = V(subnet)$name, prize = V(subnet)$prize)

edge_list <- as_edgelist(subnet)
?get.edgelist
#E(): extract edges from a network or graph object in R
edge_data <- data.frame(edge_list, weight = E(subnet)$weight)

enrichment <- res_enrichR$enrichment


write.csv(node_attributes, 'D:/RProject/Project/PCSF/CRC_PCSF_node_attributes.csv')
write.csv(edge_data, 'D:/RProject/Project/PCSF/CRC_PCSF_interactions.csv')
write.csv(enrichment, 'D:/RProject/Project/PCSF/CRC_PCSF_enrichment.csv')



