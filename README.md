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

### VI. Notes
  * Label the top directory as: td=/wsu/home/groups/piquelab/charvey/GxE
  * Scripts that call other scripts are chained together using ->
  * Sections I is performed per plate, so we will use DP1 as an example. All other plates can be processed similarly by replacing "DP1."

## Details
##### 1.) Create directories, symbolic link to fastqs, and distribute scripts for alignment
    description: create directory structure for read alignment, copy alignment/QC scripts
                 into relevant directories, & symlink to relevant .fastqs.
    script: td/derived_data/scripts/Alignment_util_makelinks.sh 
    dependencies: fastq files, named by plate number such as DP2_<things here>L2_R1.fastq 
    in: A plate number, such as DP2
    out: returns nothing

##### 2.) Align reads, merge barcodes, and perform QC 
    description:
    script: td/ Makefile 
    dependencies:
    in:
    out:

##### 3.) Read counts from all stages of filtering
    description:
    script: td/ run_all.sh -> td/ add_annotations.sh  
    dependencies:
    in:
    out:

##### 4.) Gene counts using DESeq2
    description:
    script: td/ run_all.sh -> td/ QuASAR_pipeline_masterTable.R 
    dependencies:
    in:
    out:

##### 5.) Make pileups for QuASAR
    description: 
    script: 
    in:
    out:

##### 6.) Create directories, symbolic link to pileups, and distribute scripts for QuASAR
    description: 
    script: 
    dependencies:
    in:
    out:

##### 7.) Run the QuASAR pipeline for joint genotyping and inference
    description: 
    script: td/ run_all.sh -> td/ QuASAR_pipeline_all.sh -> td/ QuASAR_pipeline.sh -> td/ QuASAR_pipeline.R
    dependencies:
    in:
    out:

##### 8.) Add gene annotations to QuASAR output
    description: 
    script: td/ add_annotations.sh 
    dependencies:
    in:
    out:

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


