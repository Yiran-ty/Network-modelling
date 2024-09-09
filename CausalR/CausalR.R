
#load CausalR
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("CausalR")

library(CausalR)


# Create Computational Causal Graph from causal network
ccg<-CreateCCG('pathway_CausalR.sif')

# Read in experimental Data -- number=Node ID from the igraph object
expData<-ReadExperimentalData('/Users/lpotarig/Library/CloudStorage/OneDrive-ImperialCollegeLondon/Student_projects/Imperial/MSc/Cancer_Informatics/Second_project/Yiran_2024/resource/CausalR_input/TF_Activity_uni.txt', ccg)
 
# Rank the predicted regulators of the experimental data from the causal network, at a max path length = 2
outrun_casualr <- runSCANR(ccg,experimentalData = expData,numberOfDeltaToScan = 4, topNumGenes = 5000)

WriteAllExplainedNodesToSifFile(outrun_casualr, ccg, expData, delta=2,
                                correctlyExplainedOnly = TRUE, quiet = TRUE)
