# GxE-pipeline

This repo contains scripts that comprise the Gene Environment Interaction (GxE) pipeline, for alignment & statistical modeling of allelic imbalances in RNA-Seq data across many cellular environments.

=================================================================

### 1. Align RNA-seq reads and perform QC
  1. align reads, merge barcodes, and perform QC
  2. read counts from all stages of filtering
  3. gene counts using DESeq2
  4. make pileups for QuASAR

### 2. QuASAR (EM) pipeline to infer individual ASE
  5.   
  6. 
### III. Process QuASAR output and prepare for MESH 
  7.)
  8.)
### IV. MESH (MCMC) to infer condition specific ASE
  9.)
  10.)
  11.)
### V. Process MESH output in preparation for multiclass logistic regression

=================================================================

#### 0.) QuASAR_util_makeLinks.sh || set up QuASAR_results_DP* directory and link to clean pileups
in:
out:
desc:

#### 1.) run_all.sh -> QuASAR_pipeline_all.sh -> QuASAR_pipeline.sh -> QuASAR_pipeline.R || run QuASAR
inp:
out: 
desc:

### 2.) run_all.sh -> add_annotations.sh || add-annotations
in:
out:
desc:

### 3.) run_all.sh -> QuASAR_pipeline_masterTable.R || combine data across a plate/cellLine into a master table
in:
out:
desc:

### 3.a.) add_annotations.sh 
in:
out:
desc:

### 4.) QuASAR_make_master.sh || concatenate all master tables from QuASAR (look into Master_table_betas_bfs.R | ./jointGenotyping_0var/data_MESH_allQuAlblfilt/Master_table_betas_bfs.R) 
in:
out:
desc:

### 5.) QuASAR_assignControls.R || find controls and assign them to all treats in the master table
in:
out:
desc:

### 6.) QuASAR_MESH_split_logFC.sh || split all data on logFC 
in:
out:
desc:

### 7.) MESH_run_all.sh || run MESH on everything
in:
out:
desc:

### 8.) master_MESH_outTable.sh || concatonate posterior CIs of configs and pi0	
in:
out:
desc:

### 9.) Master_table_betas_bfs.sh || concatenate all betas and BFs
in:
out:
desc:



## MESH


## Multinomial


     	  	   	   		     
## 1. Read alignment, quality filtering, & QC  
### Fastq -> bam 
For a new experiment, these scripts create the directory structure for read mapping and quality control of RNA-Seq reads.

## 2.
###
	
## 3.
###

## 4.
###

## 5.
###
