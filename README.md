# recount3_networks
Scripts required to replicate the analysis in - "Aggregation of publicly available RNA-seq experiments improves inference of consensus and tissue-specific Gene Co-expression Networks"

Descriptions of included directories
1. src: Directory containing R scripts required to
   1. 01_data_preprocessing: Pre-process the data and create individual groups for GTEx, TCGA, and SRA
   2. 02_data_correction: Use PCA to correct the data expression. Compute study-specific covariance matrices, weighted and unweighted aggregation of covariance matrices
   3. 03_network_inference: Infer networks for a range of penalization parameters using graphical lasso which uses study-derived covariance matrices
   4. 04_evaluate_data_aggregation: Compute log-likelihood of held-out data and F1-scores of canonical pathways for consensus networks. Compute the odds ratio of SNAP pathways for context-specific networks
   5. 05_evaluate_biological_properties: Determine the enrichment of specific GO terms among central nodes for consensus and context-specific networks. Compute the excess overlap of known gene sets among central network nodes

2. shellscripts: Bash scripts that invoke R scripts from src with corresponding inputs. The files are numbered in the order they will be invoked to produce the data subsequent scripts need.

3. plotting_notebooks: Directory that contains Rmarkdown files and html outputs to reproduce the main and supplementary figures in the paper

4. data: Directory which contains manual annotations, external validation data sets, processed expression matrices and covariance inputs for network reconstruction

5. results: Directory which contains networks inferred across a range of penalty parameters, along with network-derived statistics such as properties of the network degree distribution, computed held-out log-likelihood, F1 scores. 
