#!/bin/bash
#$ -cwd
#$ -pe smp 1
#$ -l mem_free=5G
#$ -l scratch=50G
#$ -l h_rt=24:00:00

scriptDir=/wynton/group/gladstone/biocore/projects/pfocr_pathway_enrichment_evaluation/PFOCRInPathwayAnalyses/
containerDir=/wynton/group/gladstone/biocore/containers
dataDir=/wynton/group/gladstone/biocore/projects/PFOCR/PFOCRInPathwayAnalyses

export SINGULARITY_BINDPATH="$containerDir,$scriptDir,$dataDir"


f=$1
singularity exec $containerDir/pathway_enrichment_pathway_databases_evaluation_latest.sif Rscript permute_enrichment_database.R $f

[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"