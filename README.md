# Demultiplexing Tools Testing With DEMUXAFY

running demultiplexing tools test using Demuxafy from Drew Navin (Powell's lab): 
https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Demuxlet.html

pre-demux imputation tutorial: 
https://github.com/powellgenomicslab/SNP_imputation_1000g_hg38/wiki/SNP-Genotype-Imputation-Using-1000G-hg38-Reference

--> use run_plink.sh and imputation.sh for step by step vignette, with some tips not included in the above imputation tutorial

# RUNNING THE WORKFLOW (CTP Version - use this)

to run the workflow follow the steps below and remember to run inside a `screen` session.

---

1. activate the dvc environment:  
```
conda activate /directflow/TumourProgressionGroupTemp/projects/dvc/envs/dvc_latest
```

---

2. move to top project dir: 
```
cd <PROJECT_DIR>
```

---

3. run pre-imputation snakemake step:

```
snakemake -s workflow/Snakefile_preimputation --configfile resources/testing/config/config.yaml --profile workflow/profiles/default
```

* if running from a (super) project directory: 
```
snakemake -s modules/demuxafy/workflow/Snakefile_preimputation --configfile config/config.yaml --profile modules/demuxafy/workflow/profiles/default
``` 

You should get something that looks like this:
```
Job stats:
job                                           count    min threads    max threads
------------------------------------------  -------  -------------  -------------
all                                               1              1              1
apt                                               1              1              1
apt_format_result                                 1              1              1
generate_cel_file_list_all                        1              1              1
generate_cel_genotype_sample_mappings_file        1              1              1
generate_sample_CEL_file_mappings_file            1              1              1
identify_bam_chr_encoding                         2              1              1
metadata_psam                                     1              1              1
otv_caller                                        1              1              1
ps_classification                                 1              1              1
ps_metrics                                        1              1              1
run_plink                                         1              1              1
sort_vcf                                          1              1              1
vcf_change_chr24_to_X                             1              1              1
vcf_chr_encoded                                   1              1              1
vcf_reheader                                      1              1              1
vcf_remove_position                               1              1              1
vcf_remove_unknowns                               1              1              1
vcf_rename_cel_sample_ids                         1              1              1
vcf_updated_dict                                  1              1              1
total                                            21              1              1
:
```

---

4. manually run **imputation** step (you need to actively respond to the prompts for ancestry and MAF processing) 

This step is based on the script here:
```
workflow/scripts/imputation.sh
``` 

HOWEVER, we recommend running this interactively, on the cmd line, in the following way. 

Set some container specific path bindings and add to your $PATH 
```
BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/pluginsarity
export PATH=/opt/bcftools-1.10.2/pluginsarity:$PATH
```

Define the project name - this will be used to resolve the paths for the results from the pre-imputations stage
```
PROJECT=<PROJECT_NAME>
```

Set the paths
TODO: There is too much redundancy in this path name (i.e. no need to have $PROJECT 3-times!). Let's clean this up.  
```
PROJECT_DIR=$PWD
IMPUTATION_OUTPUT_DIR=$PROJECT_DIR/results/demuxafy/$PROJECT/genotyping/$PROJECT/$PROJECT/imputation
```

Set the log output path (for HPC logs)
```
LOG=$PROJECT_DIR/logs/drmaa
mkdir -p $LOG
```

Set-up the snakemake config file required for the imputation step.

First, make a copy of the `config/imputation_config_template.yaml` file. This should be to a location (e.g., project config directory) that is relevant to the analysis being run.
```
cp config/imputation_config_template.yaml <PATH/TO/CONFIG/imputation.yaml>
```

Now, you need to edit the following fields in the template config file you've just created:
NOTE: Currently for this to work properly `PROJECT_NAME` must be the same as `config["project_id"]`
**inputs** (should look like this):  
```
inputs:
  singularity_image: resources/containers/SNP_imputation_1000g_hg38.sif # The singularity image that snakemake uses to run imputation
  plink: 
    dir: results/demuxafy/<PROJECT_NAME>/genotyping/<PROJECT_NAME>/<PROJECT_NAME>/plink ### path to directory that contains plink2 pgen, pvar and psam files (this is the ouptut from stage 1 - pre-imputation)
    pgen: actual name of the pgen file (usually "plink.pgen")
    pvar: actual name of the pvar file (usually "plink.pvar")
    psam: actual name of the psam file (usually "plink.psam")
  ref:
    dir: path to the directory containing the imputation references 
    vcf_dir: Subdirectory containing ref panel VCF
    fasta:   Genome fasta file
    genetic_map: Genetic map for phasing
    phasing_dir: Phasing subdirectory
    imputation_dir: Subdirectory with imputation references
```

**outputs**:
```
outputs:
  output_dir: results/demuxafy/<PROJECT_NAME>/genotyping/<PROJECT_NAME>/<PROJECT_NAME>/imputation # The path to the parent dir where you would like all outputs/results saved. This path must exist before running.
```

and set the `IMPUTATION_CONFIG` variable to be pointing at this `imputation.yaml` file.
```
IMPUTATION_CONFIG=<PATH/TO/CONFIG/imputation.yaml>
```

and finally, set the Snakefile that will be used to drive the imputation stage of the workflow 
```
IMPUTATION_SNAKEFILE=workflow/Snakefile_imputation
```

Check files and paths are as expected
```
echo $IMPUTATION_SNAKEFILE
echo $IMPUTATION_CONFIG
```

Do a snakemake dry-run 
```
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --dryrun --cores 1 --reason
```

You should get something that looks like this (the files will differ, but the rules shouldn't):
```
Job counts:
        count   jobs
        1       all
        1       check_sex
        1       common_snps
        1       final_pruning
        1       indiv_missingness
        1       pca_1000g
        1       pca_project
        1       pca_projection_assign
        1       prune_1000g
        9
:
```

OK, we can now run the (first part!) of the imputation workflow
```
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
JAVA_OPTS="-Xms256m -Xmx25g"

snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --profile workflow/profiles/imputation/
```

Once this all completes (successfully) then we can run the ancestry and gender update steps:
```
Ancestry_update=$IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_update_remove.tsv
Sex_update=$IMPUTATION_OUTPUT_DIR/pca_sex_checks/check_sex_update_remove.tsv

# add UPDATE in ancestry file at the endgrep line  of each line, except of the first row (starting with #FID)
sed -e '/#FID/! s/$/UPDATE/' -i $Ancestry_update

# add UPDATE in sex file at the end of each line, except of the first row (starting with #FID)
sed -e '/#FID/! s/$/UPDATE/' -i $Sex_update
```

Next, run the next part of the imputation step.
> NOTE/TODO: This is quite an unusual step for a workflow. It's very counterintuitive to use a dry-run to force updates to the parameters... It works, but we should think about if/how we can automate this or integrate into a single workflow.
```
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --dryrun --cores 1 --reason
```

At this point, you will be prompted to indicate whether you want MAF and Ancestry to be used, enter the following at the prompts: 
```
yes (at the ancestry prompt)
0   (at the MAF prompt)
```

This will then show something similar to the following rules to be run
```
Job counts:
        count   jobs
        1       all
        1       combine_vcfs_all
        2       combine_vcfs_ancestry
        2       crossmap
        44      eagle_prephasing
        2       filter_preimpute_vcf
        2       harmonize_hg38
        2       het
        2       het_filter
        44      minimac_imputation
        2       plink_to_vcf
        1       separate_indivs
        2       sort_bed
        44      split_by_chr
        2       subset_ancestry
        1       update_sex_ancestry
        2       vcf_fixref_hg38
        156
:
```

**WARNING!**  

Compare unique ancestries detected in `$IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_update_remove.tsv`, by running: 
```
cut -f4 $IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_update_remove.tsv | sort | uniq
``` 
with ancestries listed in `$IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_mafs.tsv`. 
If there are missing ancestries in the later file, **add them manually**. Keep MAF the same as indictaed for the other ancestries. 

For instance BrCa atlas data contains 4 distinct ancestries, however initially the `ancestry_mafs.tsv` file contained only 2 ancestries, while the `ancestry_update_remove.tsv` indictaed there were 4 ancestries in the dataset: 

view the file: 
```
less $IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_mafs.tsv
```

```
Ancestry    MAF
0    EAS   0.0
1    AMR   0.0
```
to check what unique ancestires were detected in the data we run: 
```
cut -f4 $IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_update_remove.tsv | sort | uniq
```  
Which will give something like the following
```
AMR
EAS
EUR
PCA_Assignment
SAS
```

Therefore, we now manually added the (2 in this example) lines for SAS and EAS ancestries to the file, resulting in: 
```
   Ancestry   MAF
0  EAS        0.0
1  AMR        0.0
2  EUR        0.0
3  SAS        0.0
```
NOTE: In the `demuxafy_test` datasets, this should only contain the following 2 ancestries, **so no changes required**.

After
1. adding UPDATE or the decision you would like to make, for the predictions, 
2. making sure ancestry_mafs.tsv contains all unique ancestries from the dataset (and updating file if necessary)
we need **can now run the second part of the imputation workflow**:
```
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
JAVA_OPTS="-Xms256m -Xmx25g"

snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --profile workflow/profiles/imputation/
```
This will generate a vcf file that is used as an input for the demux tools in the post-imputation steps (see below).

----

5. run the post-imputation part of the demuxafy snakemake: 
```
TODO: REMOVE: Not needed as we can use the same snakemake env throughout
#conda deactivate
#conda activate /directflow/TumourProgressionGroupTemp/projects/dvc/envs/dvc_stable

snakemake -s workflow/Snakefile_postimputation --configfile resources/testing/config/config.yaml --profile workflow/profiles/default
```

This will result in something like the following (on test-dataset)
```
Job stats:
job                                        count    min threads    max threads
---------------------------------------  -------  -------------  -------------
add_vireosnp_to_combined_results               2              1              1
all                                            1              1              1
combine_cell_assignments                       1              1              1
combine_cell_stats                             1              1              1
combine_singlet_stats                          1              1              1
demuxafy                                       2              1              1
demuxafy_count_cells_AtLeastHalfSinglet        2              1              1
demuxlet                                       2              1              1
dsc_pileup                                     2              1              1
dsc_pileup_demuxlet                            2              1              1
freemuxlet                                     2              1              1
freemuxlet_summary                             2              1              1
pdf_to_png                                     2              1              1
report_demuxafy                                1              1              1
souporcell                                     2              1              1
tabix_vcf                                      2              1              1
vcf_filter                                     1              1              1
vcf_filter2                                    1              1              1
vireo                                          2              1              1
vireosnp                                       2              1              1
wrangle_demuxafy_output                        2              1              1
zip_vcf                                        2              1              1
total                                         37              1              1
:
```

---
# STRUCTURE OF THE WORKFLOW

As described above, the workflow is currently split into 3 stages:

- pre-imputation
- imputation
- post-imputation

The rule DAG for each of these stages is as follows:  
NOTE: generated using commands like this:
```
snakemake -s workflow/Snakefile_postimputation --configfile resources/testing/config/config_renamed.yaml --profile conda_profile  --forceall --rulegraph | dot -Tpng > temp/dag.imputation-stage_3.png
```

## pre-imputation

![pre-imputation](resources/images/dag.pre_imputation-stage_1.png)

## imputation

![imputation](resources/images/dag.imputation-stage_2.png)

## post-imputation

![post-imputation](resources/images/dag.post_imputation-stage_3.png)

------

# RUN THE TEST DATASET

1a. unprotect some files that cause errors if they are sym-links
```
dvc unprotect resources/reference/hg38/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
dvc unprotect resources/genomes/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa.fai
```

1. run the pre-imputation step:

```
snakemake -s workflow/Snakefile_preimputation --configfile resources/testing/config/config.mixed_pool12.missing_PBMCs.yaml --profile workflow/profiles/default
```

You should see this (use a `--dry-run``):
```
Job stats:
job                                           count    min threads    max threads
------------------------------------------  -------  -------------  -------------
all                                               1              1              1
apt                                               1              1              1
apt_format_result                                 1              1              1
generate_cel_file_list_all                        1              1              1
generate_cel_genotype_sample_mappings_file        1              1              1
generate_sample_CEL_file_mappings_file            1              1              1
identify_bam_chr_encoding                         2              1              1
metadata_psam                                     1              1              1
otv_caller                                        1              1              1
ps_classification                                 1              1              1
ps_metrics                                        1              1              1
run_plink                                         1              1              1
sort_vcf                                          1              1              1
vcf_change_chr24_to_X                             1              1              1
vcf_chr_encoded                                   1              1              1
vcf_reheader                                      1              1              1
vcf_remove_position                               1              1              1
vcf_remove_unknowns                               1              1              1
vcf_rename_cel_sample_ids                         1              1              1
vcf_updated_dict                                  1              1              1
total                                            21              1              1
:
```

---


2. manually run **imputation** step (you need to actively respond to the prompts for ancestry and MAF processing) 

Set some container specific path bindings and add to your $PATH 
```
BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/pluginsarity
export PATH=/opt/bcftools-1.10.2/pluginsarity:$PATH
```

Define the project name - this will be used to resolve the paths for the results from the pre-imputations stage
```
PROJECT=test_mixed_pool12.missing_PBMCs
```

Set the paths 
```
PROJECT_DIR=$PWD
IMPUTATION_OUTPUT_DIR=$PROJECT_DIR/results/demuxafy/$PROJECT/genotyping/$PROJECT/$PROJECT/imputation
```

Set the log output path (for HPC logs)
```
LOG=$PROJECT_DIR/logs/drmaa
mkdir -p $LOG
```

Set-up the snakemake config file required for the imputation step
and set the `IMPUTATION_CONFIG` variable to be pointing at this `imputation.yaml` file.
```
IMPUTATION_CONFIG=resources/testing/config/imputation.mixed_pool12.missing_PBMCs.yaml
```

and finally, set the Snakefile that will be used to drive the imputation stage of the workflow 
```
IMPUTATION_SNAKEFILE=workflow/Snakefile_imputation
```

Check files and paths are as expected
```
echo $IMPUTATION_SNAKEFILE
echo $IMPUTATION_CONFIG
```

Do a snakemake dry-run 
```
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --dryrun --cores 1 --reason
```

You should get something that looks like this:
```
Using these files from the plink directory results/demuxafy/test_mixed_pool12.missing_PBMCs/genotyping/test_mixed_pool12.missing_PBMCs/test_mixed_pool12.missing_PBMCs/plink as input:
The pgen file: results/demuxafy/test_mixed_pool12.missing_PBMCs/genotyping/test_mixed_pool12.missing_PBMCs/test_mixed_pool12.missing_PBMCs/plink/plink.pgen
The pvar file: results/demuxafy/test_mixed_pool12.missing_PBMCs/genotyping/test_mixed_pool12.missing_PBMCs/test_mixed_pool12.missing_PBMCs/plink/plink.pvar
The psam file: results/demuxafy/test_mixed_pool12.missing_PBMCs/genotyping/test_mixed_pool12.missing_PBMCs/test_mixed_pool12.missing_PBMCs/plink/plink.psam

Found all the required refernce files in resources/reference/hg38

Building DAG of jobs...
Job stats:
job                      count    min threads    max threads
---------------------  -------  -------------  -------------
all                          1              1              1
check_sex                    1              1              1
common_snps                  1              1              1
final_pruning                1              1              1
indiv_missingness            1              1              1
pca_1000g                    1              1              1
pca_project                  1              1              1
pca_projection_assign        1              1              1
prune_1000g                  1              1              1
total                        9              1              1

:
```

OK, we can now run the (first part!) of the imputation workflow
```
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
JAVA_OPTS="-Xms256m -Xmx25g"

snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --profile workflow/profiles/imputation/
```

Once this all completes (successfully) then we can run the ancestry and gender update steps:
```
Ancestry_update=$IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_update_remove.tsv
Sex_update=$IMPUTATION_OUTPUT_DIR/pca_sex_checks/check_sex_update_remove.tsv

# add UPDATE in ancestry file at the endgrep line  of each line, except of the first row (starting with #FID)
sed -e '/#FID/! s/$/UPDATE/' -i $Ancestry_update

# add UPDATE in sex file at the end of each line, except of the first row (starting with #FID)
sed -e '/#FID/! s/$/UPDATE/' -i $Sex_update
```

Next, run the next part of the imputation step.
```
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --dryrun --cores 1 --reason
```

At this point, you will be prompted to indicate whether you want MAF and Ancestry to be used, enter the following at the prompts: 
```
yes (at the ancestry prompt)
0   (at the MAF prompt)
```

This will then show something similar to the following rules to be run
```
Job counts:
        count   jobs
        1       all
        1       combine_vcfs_all
        2       combine_vcfs_ancestry
        2       crossmap
        44      eagle_prephasing
        2       filter_preimpute_vcf
        2       harmonize_hg38
        2       het
        2       het_filter
        44      minimac_imputation
        2       plink_to_vcf
        1       separate_indivs
        2       sort_bed
        44      split_by_chr
        2       subset_ancestry
        1       update_sex_ancestry
        2       vcf_fixref_hg38
        156
:
```

**WARNING!**  

to check what unique ancestires were detected in the data we run: 
```
cut -f4 $IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_update_remove.tsv | sort | uniq
```  
Which should give the following
```
EAS
EUR
PCA_Assignment
```

NOTE: In the `demuxafy_test` datasets, this should only contain these 2 ancestries, **so no changes required**.

we **can now run the second part of the imputation workflow**:
```
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
JAVA_OPTS="-Xms256m -Xmx25g"

snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --profile workflow/profiles/imputation/
```
This will generate a vcf file that is used as an input for the demux tools in the post-imputation steps (see below).

----

TODO: from here

3. run the post-imputation part of the demuxafy snakemake: 

```
snakemake -s workflow/Snakefile_postimputation --configfile resources/testing/config/config.mixed_pool12.missing_PBMCs.yaml --profile workflow/profiles/default
```

This will result in the following jobs to run:
```
Job stats:
job                                        count    min threads    max threads
---------------------------------------  -------  -------------  -------------
add_vireosnp_to_combined_results               2              1              1
all                                            1              1              1
combine_cell_assignments                       1              1              1
combine_cell_stats                             1              1              1
combine_singlet_stats                          1              1              1
demuxafy                                       2              1              1
demuxafy_count_cells_AtLeastHalfSinglet        2              1              1
demuxlet                                       2              1              1
demuxlet_summary                               2              1              1
dsc_pileup_demuxlet                            2              1              1
dsc_pileup_freemuxlet                          2              1              1
freemuxlet                                     2              1              1
freemuxlet_summary                             2              1              1
pdf_to_png                                     2              1              1
report_demuxafy                                1              1              1
souporcell                                     2              1              1
tabix_vcf                                      2              1              1
vcf_filter                                     1              1              1
vcf_filter2                                    1              1              1
vireo                                          2              1              1
vireosnp                                       2              1              1
wrangle_demuxafy_output                        2              1              1
zip_vcf                                        2              1              1
total                                         39              1              1

:
```

If doing any development work, run the test compare script to see if the test has been successful
```
./resources/testing/run_test_compare.mixed_pool12.missing_PBMCs.sh > resources/testing/run_test_compare.mixed_pool12.missing_PBMCs.output.log 2>&1
```
and check the results output `log`` file to check all steps PASSED. 