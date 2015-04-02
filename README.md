# GxE-pipeline

This repo contains scripts that comprise the Gene Environment Interaction (GxE) pipeline, for alignment & statistical modeling of allelic imbalances in RNA-Seq data across many cellular environments.

=================================================================
## I. Align RNA-seq reads and perform QC
    0. Create directories, symbolic link to fastqs, and distribute scripts for alignment
    1. Align reads, merge barcodes, and perform QC
    2. Read counts from all stages of filtering
    3. Gene counts using DESeq2
    4. Make pileups for QuASAR
## II. QuASAR (EM) pipeline to infer individual ASE
    5. Create directories, symbolic link to pileups, and distribute scripts for QuASAR
    6. Run the QuASAR pipeline for joint genotyping and inference   
    7. Add gene annotations to QuASAR output
## III. Process QuASAR output and prepare for MESH 
    8. Combine all QuASAR output into a master table
    9. Add logFC annotations (Greg's script)
    10. Create directories and distribute scripts for MESH	
    11. Split data for MESH and distribute to analysis directory
## IV. MESH (MCMC) to infer condition specific ASE
    9. Run the MESH pipeline
    10. Calculate posteriors and bayes factors
    11. Plot condition specific ASE
    12. Prep data for Multinomial
## V. Process MESH output in preparation for multiclass logistic regression
    13. Run multiclass logistic regression
    14. Plot regression results

=================================================================
##### 0.) Alignment_util_makelinks.sh || set up QuASAR_results_DP* directory and link to clean pileups
    script: 
    description: 
    dependencies:
    in:
    out:

##### 1.) run_all.sh -> QuASAR_pipeline_all.sh -> QuASAR_pipeline.sh -> QuASAR_pipeline.R || run QuASAR
    script: 
    description: 
    dependencies:
    in:
    out:

##### 2.) run_all.sh -> add_annotations.sh || add-annotations
    script: 
    description: 
    dependencies:
    in:
    out:

##### 3.) run_all.sh -> QuASAR_pipeline_masterTable.R || combine data across a plate/cellLine into a master table
    script: 
    description: 
    dependencies:
    in:
    out:

##### 3.a.) add_annotations.sh 
    script: 
    description: 
    dependencies:
    in:
    out:

##### 4.) QuASAR_make_master.sh || concatenate all master tables from QuASAR (look into Master_table_betas_bfs.R | ./jointGenotyping_0var/data_MESH_allQuAlblfilt/Master_table_betas_bfs.R) 
    script: 
    description: 
    dependencies:
    in:
    out:

##### 5.) QuASAR_assignControls.R || find controls and assign them to all treats in the master table
    script: 
    description: 
    dependencies:
    in:
    out:

##### 6.) QuASAR_MESH_split_logFC.sh || split all data on logFC 
    script: 
    description: 
    dependencies:
    in:
    out:

##### 7.) MESH_run_all.sh || run MESH on everything
    script: 
    description: 
    dependencies:
    in:
    out:

##### 8.) master_MESH_outTable.sh || concatonate posterior CIs of configs and pi0	
    script: 
    description: 
    dependencies:
    in:
    out:


##### 9.) Master_table_betas_bfs.sh || concatenate all betas and BFs
    script: 
    description: 
    dependencies:
    in:
    out:
