# =============
# IMPUTATION
# ============= 

# activate snakemake env
conda activate wg1_snakemake
# Container specific binding
BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/pluginsarity
#add BCFTOOLS_PLUGINS to your $PATH 
export PATH=/opt/bcftools-1.10.2/pluginsarity:$PATH

PROJECT=BrCa_Pilot2
# TODO: [IMPORTANT] make configurable outside script
#PROJECT_DIR=/share/ScratchGeneral/danrod/workflows/demuxafy/demultiplexing_tools_testing_with_DEMUXAFY
PROJECT_DIR=$PWD
IMPUTATION_OUTPUT_DIR=$PROJECT_DIR/results/demuxafy/genotyping/$PROJECT/$PROJECT/imputation
# NOTE/TODO: what is logged here?
LOG=$IMPUTATION_OUTPUT_DIR
mkdir -p $LOG

# NOTE: could probably re-name this to be the config.yaml for the snakemake pipeline
# IMPORTANT: at the moment the paths to bam file, barcode, vcf and plink are all given in the PreImputation.yaml file! 
# TO DO: move the input variables from yaml file
# substitute variables

# TODO: streamline this - it seems unnecessary to have a sepoarate step to edit a config like this
IMPUTATION_CONFIG=image/PreImputation_edited.yaml
# IMPUTATION_CONFIG_original=image/PreImputation.yaml
# IMPUTATION_CONFIG=image/PreImputation_edited.yaml
# sed 's|PROJECT_DIR|'"$PROJECT_DIR"'|g' $IMPUTATION_CONFIG_original > $IMPUTATION_CONFIG
# sed -i 's|PROJECT_NAME|'"$PROJECT"'|g' $IMPUTATION_CONFIG 
# sed -i 's|OUTPUT_DIR|'"$LOG"'|g' $IMPUTATION_CONFIG

#IMPUTATION_CONFIG=/directflow/TumourProgressionGroupTemp/projects/tmp/beakie/brca_atlas/brca_atlas-data/modules/demuxafy/config/demuxafy/PreImputation.yaml
#IMPUTATION_SNAKEFILE=/directflow/TumourProgressionGroupTemp/projects/tmp/beakie/demuxafy/image/Snakefile
IMPUTATION_SNAKEFILE=image/Snakefile

# LOG dir specified above, needs to be the same as the output_dir in the imputation_config 

####CHECK YOUR FILES AND PATH IS CORRECT####
echo $IMPUTATION_SNAKEFILE
echo $IMPUTATION_CONFIG

# run dry-run 
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --dryrun --cores 1 --reason

# run!!!
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
JAVA_OPTS="-Xms256m -Xmx25g"
#export BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/plugins

# TODO: run via drmaa (actually abstract this from here and make system agnostic)
#nohup snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --rerun-incomplete --jobs 20 --use-singularity --restart-times 2 --keep-going --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -P TumourProgression -e $LOG -o $LOG -j y -V" > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --rerun-incomplete --jobs 20 --use-singularity --restart-times 2 --keep-going --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -P TumourProgression -e $LOG -o $LOG -j y -V" > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log

#manually chnaged #UPDATE to /directflow/TumourProgressionGroupTemp/projects/tmp/beakie/demuxafy/results/combined/pca_sex_checks files 
# after initial run add UPFATE to initial column

Ancestry_update=$IMPUTATION_OUTPUT_DIR/pca_sex_checks/ancestry_update_remove.tsv
Sex_update=$IMPUTATION_OUTPUT_DIR/pca_sex_checks/check_sex_update_remove.tsv

# add UPDATE in ancestry file at the endgrep line  of each line, except of the first row (starting with #FID)
sed -e '/#FID/! s/$/UPDATE/' -i $Ancestry_update
# this was an attempt to add extra tab (\t), but the file is already a tsv, so it adds it automatically
#sed -e '/#FID/! s/$/\tUPDATE/' -i test.tsv

# add UPDATE in sex file at the end of each line, except of the first row (starting with #FID)
sed -e '/#FID/! s/$/UPDATE/' -i $Sex_update

# after adding UPDATE or the decision you would like to make for the predictions run the pipleine again 
# NOTE/TODO: This is a really strange step to the workflow. It's very counterintuitive to use a dry-run to force updates to the parameters... Think about if/how we can automate this or integrate into a single workflow..
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --dryrun --cores 1 --reason

# here you need to indictae whether you want MAF and Ancestry to be used 
# do yes for ancestry and 0 for MAF
# so in the dry-run it asks you your preference for ancestry imputation and MAF

# run!!!
export OMP_NUM_THREADS=1
export USE_SIMPLE_THREADED_LEVEL3=1
JAVA_OPTS="-Xms256m -Xmx25g"
#export BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/plugins

# TODO: add this to a drmaa based submission profile
snakemake --snakefile $IMPUTATION_SNAKEFILE --configfile $IMPUTATION_CONFIG --rerun-incomplete --jobs 20 --use-singularity --restart-times 2 --keep-going --cluster "qsub -S /bin/bash -q short.q -r yes -pe smp {threads} -l tmp_requested={resources.disk_per_thread_gb}G -l mem_requested={resources.mem_per_thread_gb}G -P TumourProgression -e $LOG -o $LOG -j y -V" > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log

# nohup \
#     snakemake \
#       --snakefile $IMPUTATION_SNAKEFILE \
#       --configfile $IMPUTATION_CONFIG \
#       --rerun-incomplete \
#       --jobs 20 \
#       --use-singularity \
#       --restart-times 2 \
#       --keep-going \
#       --cluster \
#           "qsub -S /bin/bash \
#           -q short.q \
#           -r yes \
#           -pe smp {threads} \
#           -l tmp_requested={resources.disk_per_thread_gb}G \
#           -l mem_requested={resources.mem_per_thread_gb}G \
#           -P TumourProgression \
#           -e $LOG \
#           -o $LOG \
#           -j y \
#           -V" \
#     > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &


# this will generate a vcf file that can be used in as an input for demux tools after filtering: 

# NOTE (DLR)/TODO: Is this a comment related to a specific analysis or should it always be run incase some CEL files didn't make it?
# not all CEL files made it. 
# re-run to include individuals that didnt pass the threshold

SIF=$PROJECT_DIR/image/SNP_imputation_1000g_hg38.sif
singularity exec --bind /share,/directflow $SIF plink2 --missing --pfile $PROJECT_DIR/results/demuxafy/genotyping/BrCa_Pilot2/BrCa_Pilot2/plink/plink > $PROJECT_DIR/results/demuxafy/genotyping/BrCa_Pilot2/BrCa_Pilot2/plink/test

# =============
# RUN postimputation.smk
# ============= 