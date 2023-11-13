# recount3_networks
Scripts required to replicate the analysis in - "Aggregation of publicly available RNA-seq experiments improves inference of consensus and tissue-specific Gene Co-expression Networks"

Note: Prior to running the scripts, please ensure that the homeDir variable is correctly populated with the location of the directory that you have cloned this repo into. 

Descriptions of included directories
1. src: Directory containing R scripts required to
   1. 01_data_preprocessing: Pre-process the data and create individual groups for GTEx, TCGA, and SRA
   2. 02_data_correction: Use PCA to correct the data expression. Compute study-specific covariance matrices, weighted and unweighted aggregation of covariance matrices
   3. 03_network_inference: Infer networks for a range of penalization parameters using graphical lasso which uses study-derived covariance matrices
   4. 04_evaluate_data_aggregation: Compute log-likelihood of held-out data and F1-scores of canonical pathways for consensus networks. Compute the odds ratio of SNAP pathways for context-specific networks
   5. 05_evaluate_biological_properties: Determine the enrichment of specific GO terms among central nodes for consensus and context-specific networks. Compute the excess overlap of known gene sets among central network nodes

2. shellscripts: Bash scripts that invoke R scripts from src with corresponding inputs. The files are numbered in the order they will be invoked to produce the data subsequent scripts need.

3. plotting_notebooks: Directory that contains Rmarkdown files and html outputs to reproduce the main and supplementary figures in the paper.
   1. Data_overview.rmd: Used to replicate Main Figure 1, and supplementary figures 1 and 2, describing the properties of the pre-processed and aggregated RNA-seq data 
   2. Modeling_confounder.rmd: Used to understand the major axes of variation in observed gene expression as a function of varying levels of data aggregation. Used to replicate supplementary figures 3 and 4.
   3. Evaluate_data_aggregation.rmd: Used to determine the improvement in held-out log-likelihood and F1 scores of networks with data aggregation in GTEx and SRA splits. Generates Main Fig. 2, supplementary figures 5 - 9
   4. Evaluate_data_aggregation_consensus_context_specific.rmd: Determine the improvement in the inference of consensus and context-specific networks belonging to blood and CNS. Generates main Fig. 3. 
   5. Enrichment_GO_terms.rmd: Visualize the enrichment of GO terms with an increase in network centrality in consensus and context-specific networks. Main Figures 4 and 5
   6. Excess_overlap_gene_sets.rmd: Visualize the excess overlap of functionally defined and evolutionarily conserved gene sets among network nodes binned by closeness centrality Main Figures 4 and 5
   7. Properties_of_context_specific_networks.rmd: Visualizing enrichment of tissue-specific and general transcription factors among network central nodes. Visualize matrix factorization of context-specific networks based on shared hubs and edges.
   8. Heritability_enrichment_consensus.rmd: Visualize the heritability enrichment and coefficient of centrality annotations derived from universal and non-cancerous consensus networks by meta-analyzing either across 42 independent traits or 219 UKBB traits.
   9. Heritability_enrichment_context_specific.rmd: Visualize the heritability enrichment and coefficient of centrality annotations derived from blood and CNS networks by meta-analyzing either across concordant traits for each context. 

5. data: Directory which contains manual annotations, external validation data sets, processed expression matrices and covariance inputs for network reconstruction

6. results: Directory which contains networks inferred across a range of penalty parameters, along with network-derived statistics such as properties of the network degree distribution, computed held-out log-likelihood, F1 scores. 
