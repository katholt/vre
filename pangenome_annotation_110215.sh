#!/bin/bash
#SBATCH -p main
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16384
#SBATCH --time=05:00:00

module load prokka/1.10

prokka --outdir EFAE_pan_all_110215 --cpus 1 --locustag EFAE --genus Enterococcus --species faecium --usegenus --proteins /vlsci/VR0082/shared/resistanceDB_srst2/resistance.fasta /vlsci/VR0082/shared/data/enterococcus/faecium/refgenomes/pangenome/VRE_pangenome_renamed.fasta

jhawkey@barcoo /vlsci/VR0082/shared/data/enterococcus/faecium/refgenomes/pangenome $ sbatch prokka_pan_all.txt
Submitted batch job 3157105

#!/bin/bash
#SBATCH -p main
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=16384
#SBATCH --time=05:00:00

module load prokka/1.10

prokka --cpus 1 --locustag EFAE --genus Enterococcus --species faecium --usegenus --proteins /vlsci/VR0082/shared/resistanceDB_srst2/resistance.fasta /vlsci/VR0082/shared/data/enterococcus/faecium/refgenomes/pangenome/VRE_pangenome_part_renamed.fasta

jhawkey@barcoo /vlsci/VR0082/shared/data/enterococcus/faecium/refgenomes/pangenome $ sbatch prokka_pan_part.txt
Submitted batch job 3157106

