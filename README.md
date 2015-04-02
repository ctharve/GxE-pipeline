# GxE-pipeline

This repo contains scripts that comprise the Gene Environment Interaction (GxE) pipeline, for alignment & statistical modeling of allelic imbalances in RNA-Seq data across many cellular environments.

## I. Align RNA-seq reads and perform QC
    1. Create directories, symbolic link to fastqs, and distribute scripts for alignment
    2. Align reads, merge barcodes, and perform QC
    3. Read counts from all stages of filtering
    4. Gene counts using DESeq2
    5. Make pileups for QuASAR
## II. QuASAR (EM) pipeline to infer individual ASE
    6. Create directories, symbolic link to pileups, and distribute scripts for QuASAR
    7. Run the QuASAR pipeline for joint genotyping and inference   
    8. Add gene annotations to QuASAR output
## III. Process QuASAR output and prepare for MESH 
    9. Combine all QuASAR output into a master table
    10. Add logFC annotations (Greg's script)
    11. Create directories and distribute scripts for MESH	
    12. Split data for MESH and distribute to analysis directory
## IV. MESH (MCMC) to infer condition specific ASE
    13. Run the MESH pipeline
    14. Calculate posteriors and bayes factors
    15. Plot condition specific ASE
    16. Prep data for Multinomial
## V. Process MESH output in preparation for multiclass logistic regression
    17. Run multiclass logistic regression
    18. Plot regression results
## VI. Notes
    * Label the top directory as: td=/wsu/home/groups/piquelab/charvey/GxE
    * Scripts that call other scripts are chained together using ->
=================================================================
##### 1. Create directories, symbolic link to fastqs, and distribute scripts for alignment
    script: td/derived_data/scripts/Alignment_util_makelinks.sh 
    description: 
    dependencies:
    in:
    out:

##### 2. Align reads, merge barcodes, and perform QC 
    script: 
    description: 
    dependencies:
    in:
    out:

##### 3. Read counts from all stages of filtering
    script: td/ run_all.sh -> td/ add_annotations.sh 
    description: 
    dependencies:
    in:
    out:

##### 4. Gene counts using DESeq2
    script: td/ run_all.sh -> td/ QuASAR_pipeline_masterTable.R
    description: 
    dependencies:
    in:
    out:

##### 5. Make pileups for QuASAR
    script: 
    description: 
    dependencies:
    in:
    out:

##### 6. Create directories, symbolic link to pileups, and distribute scripts for QuASAR
    script: 
    description: 
    dependencies:
    in:
    out:

##### 7. Run the QuASAR pipeline for joint genotyping and inference
    script: td/ run_all.sh -> td/ QuASAR_pipeline_all.sh -> td/ QuASAR_pipeline.sh -> td/ QuASAR_pipeline.R
    description: 
    dependencies:
    in:
    out:

##### 8. Add gene annotations to QuASAR output
    script: td/ add_annotations.sh 
    description: 
    dependencies:
    in:
    out:

##### 9. Combine all QuASAR output into a master table
    script: QuASAR_make_master.sh || concatenate all master tables from QuASAR (look into Master_table_betas_bfs.R | ./jointGenotyping_0var/data_MESH_allQuAlblfilt/Master_table_betas_bfs.R) 
    description: 
    dependencies:
    in:
    out:

##### 10. Add logFC annotations (Greg's script)
    script: 
    description: 
    dependencies:
    in:
    out:


##### 11. Create directories and distribute scripts for MESH
    script: 
    description: 
    dependencies:
    in:
    out:

##### 12. Split data for MESH and distribute to analysis directory
    script: QuASAR_assignControls.R || find controls and assign them to all treats in the master table
    script1: QuASAR_MESH_split_logFC.sh || split all data on logFC 
    description: 
    dependencies:
    in:
    out:

##### 13. Run the MESH pipeline
    script: td/ MESH_run_all.sh || run MESH on everything
    description: 
    dependencies:
    in:
    out:

##### 14. Calculate posteriors and bayes factors
    script: td/ Master_table_betas_bfs.sh || concatenate all betas and BFs 
    script1: master_MESH_outTable.sh || concatonate posterior CIs of configs and pi0
    description: 
    dependencies:
    in:
    out:

##### 15. Plot condition specific ASE
    script: 
    description: 
    dependencies:
    in:
    out:

##### 16. Prep data for Multinomial
    script: 
    description: 
    dependencies:
    in:
    out:

##### 17. Run multiclass logistic regression
    script: 
    description: 
    dependencies:
    in:
    out:

##### 18. Plot regression results 
    script: 
    description: 
    dependencies:
    in:
    out:


