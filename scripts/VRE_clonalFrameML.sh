#!/bin/bash

#SBATCH --job-name=CFML_VRE

#SBATCH --time=48:00:00

#SBATCH --ntasks=1

#SBATCH --mem-per-cpu=32768



module load clonalframeml-gcc/1.8

ClonalFrameML 95_gene_cons_RAxML_bipartitions.rooted.tree common_95_variant_OG_removed_CP006620_wg.aln emsim100.output -emsim 100 > VRE_emsim100_0704.log 