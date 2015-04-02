# GxE-pipeline
This repo contains scripts that comprise the Gene Environment Interaction (GxE) pipeline, for alignment & statistical modeling of allelic imbalances in RNA-Seq data across many cellular environments.

## Outline
### I. Align RNA-seq reads and perform QC
  * 1.) Create directories, symbolic link to fastqs, and distribute scripts for alignment
  * 2.) Align reads, merge barcodes, and perform QC
  * 3.) Read counts from all stages of filtering
  * 4.) Gene counts using DESeq2
  * 5.) Make pileups for QuASAR

### II. QuASAR (EM) pipeline to infer individual ASE
  * 6.) Create directories, symbolic link to pileups, and distribute scripts for QuASAR
  * 7.) Run the QuASAR pipeline for joint genotyping and inference   
  * 8.) Add gene annotations to QuASAR output

### III. Process QuASAR output and prepare for MESH 
  * 9.) Combine all QuASAR output into a master table
  * 10.) Add logFC annotations (Greg's script)
  * 11.) Create directories and distribute scripts for MESH	
  * 12.) Split data for MESH and distribute to analysis directory

### IV. MESH (MCMC) to infer condition specific ASE
  * 13.) Run the MESH pipeline
  * 14.) Calculate posteriors and bayes factors
  * 15.) Plot condition specific ASE
  * 16.) Prep data for Multinomial

### V. Process MESH output in preparation for multiclass logistic regression
  * 17.) Run multiclass logistic regression
  * 18.) Plot regression results

## Notes
  * Label the top directory as: td=/wsu/home/groups/piquelab/charvey/GxE
  * Scripts that call other scripts are chained together using ->
  * I is performed per plate, so we will use DP1 as an example. All other plates can be processed similarly by replacing DP1 with the desired plate.

## Details
##### 1.) Create directories, symbolic link to fastqs, and distribute scripts for alignment
    description: create directory structure for read alignment, copy alignment/QC scripts
                 into relevant directories, & symlink to relevant .fastqs.
    script: td/derived_data/scripts/Alignment_util_makelinks.sh 
    dependencies: fastq files, named by plate number such as DP1_<things here>L1_R1.fastq 
    in: A plate number, such as DP1
    out: null

##### 2.) Align reads, merge barcodes, and perform QC 
    description: align fastq files, merge barcodes, quality filter, then remove duplicates. Do read counts at each stage of filtering
    script: td/derived_data/DP1/bams/Makefile 
    dependencies: the folder ../fastqs with symlinked fastqs.
    in: null
    out: clean bam files & .txt files with read counts from each stage of filtering

##### 3.) Read counts from all stages of filtering
    description: creates a histogram & table from read counts at each stage of .bam processing
    script: td/derived_data/DP1/counts/QC counts_QC_logs.sh -> counts_QC_logs.R  
    dependencies: 2.) is succesfully completed
    in: A plate number, such as DP1
    out: DP1_QC_counts.pdf & DP1_QC_counts.txt

##### 4.) Gene counts using DESeq2
    description: use DESeq2 to measure gene expression, at the transcript level, from processed bams
    script: td/derived_data/DP1/counts/GC counts_DEG.sh -> counts_DEG.R
    dependencies: 2.) is succesfully completed
    in: A plate number, such as DP1
    out: DP1.data.gz & DP1.data.Rd

##### 5.) Make pileups for QuASAR
    description: make pileup/bed files from each clean bam file
    script: td/derived_data/DP1/pileups/Makefile1
    dependencies: 2.) is succesfully completed
    in: null
    out: clean DP1*.pileup.clean.bed.gz files 

##### 6.) Create directories, symbolic link to pileups, and distribute scripts for QuASAR
    description: create directory structure for QuASAR pipeline, copy QuASAR/MESH scripts
                 into relevant directories, & symlink to relevant *.clean.bed.gz files.
    script: td/jointGenotyping/scripts/QuASAR_util_makeLinks.sh
    dependencies: 5.) is succesfully completed
    in: A plate number, such as DP1
    out: null

##### 7.) Run the QuASAR pipeline for joint genotyping and inference
    description: a full QuASAR analysis for each plate/cellLine
    script: td/jointGenotyping/QuASAR_results_DP1/ QuASAR_pipeline_all.sh ->  QuASAR_pipeline.R
    dependencies: 6.) is succesfully completed & a properly formatted covariate file in td/derived_data/covariates 
    in: A plate number, such as DP1, & the name to the analysis script, QuASAR_pipeline.R
    out: QuASAR output for each plate:cellLine:treatment, plate:cellLine, QQ plots, and manhattan plots

##### 8.) Add gene annotations to QuASAR output
    description: add annotations of each gene to every heterozygote SNP identified with QuASAR
    script: td/jointGenotyping/QuASAR_results_DP1/output/add_annotations.sh 
    dependencies: 7.) is succesfully completed
    in: null
    out: QuASAR output tables with gene name annotation on each SNP

##### 9. Combine all QuASAR output into a master table
    description: 
    script: QuASAR_make_master.sh || concatenate all master tables from QuASAR (look into Master_table_betas_bfs.R | ./jointGenotyping_0var/data_MESH_allQuAlblfilt/Master_table_betas_bfs.R) 
    dependencies:
    in:
    out:

##### 10. Add logFC annotations (Greg's script)
    description:  
    script: 
    dependencies:
    in:
    out:


##### 11. Create directories and distribute scripts for MESH
    description: 
    script: 
    dependencies:
    in:
    out:

##### 12. Split data for MESH and distribute to analysis directory
    description: 
    script: QuASAR_assignControls.R || find controls and assign them to all treats in the master table
    script1: QuASAR_MESH_split_logFC.sh || split all data on logFC 
    dependencies:
    in:
    out:

##### 13. Run the MESH pipeline
    description: 
    script: td/ MESH_run_all.sh || run MESH on everything
    dependencies:
    in:
    out:

##### 14. Calculate posteriors and bayes factors
    description: 
    script: td/ Master_table_betas_bfs.sh || concatenate all betas and BFs 
    script1: master_MESH_outTable.sh || concatonate posterior CIs of configs and pi0
    dependencies:
    in:
    out:

##### 15. Plot condition specific ASE
    description: 
    script: 
    dependencies:
    in:
    out:

##### 16. Prep data for Multinomial
    description: 
    script: 
    dependencies:
    in:
    out:

##### 17. Run multiclass logistic regression
    description: 
    script: 
    dependencies:
    in:
    out:

##### 18. Plot regression results 
    description: 
    script: 
    dependencies:
    in:
    out:


