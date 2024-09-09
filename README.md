# Network-modelling and Structural comparison
### 1.HMI
DMI_terminal: the MicrobioLink pipeline for Host-Microbe interaction prediction

IUPred: the IUPred tool to only retains the motifs within intrinsically disordered regions
### 2. bacteria_domain_structure
download_bacterial_proteins: download the bacteria proteins from Uniprot Proteome database
### 3. enrichment_analysis
enrichr_id_database_ranking: the functional enrichment analysis for TieDIE and CARNIVAL outcome using Reactome database.
### 4. location_filter_human
get_human_fasta: filter the human protein based on the location, download the human protein fasta from UniProt database.
### 5. tiedie
The TieDIE pipeline

TieDie_input_process: data processing to prepare the input for TieDIE

TieDie: TieDIE main algorithm

TieDie_output_process: further processing of the tieDIE output, connecting the intermediate signalling network with upstream bacteria proteins and downstream DEGs.
### 6. CARNIVAL
carnival.R: CARNIVAL algorithm

carnival_output_processing: further processing of the CARNIVAL output, connecting the intermediate signalling network with upstream bacteria proteins and downstream DEGs.
### 7. PCSF
PCSF_CRC.R: PCSF analysis and the enrichment analysis for the PCSF output

### 8. CausalR

