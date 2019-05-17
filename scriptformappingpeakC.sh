#!/bin/bash
#SBATCH --job-name=PeakC131
#SBATCH --time=11:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=40gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

inputfolder=/groups/umcg-gdio/tmp04/umcg-ihidding/4C/peakC/inputfolder

module load BWA
module load BEDTools
module load PerlPlus
export PERL5LIB=/home/umcg-ihidding/perl5/lib/perl5/
numb=0
mkdir /groups/umcg-gdio/tmp04/umcg-ihidding/4C/peakC/storedoutput3

for f in $inputfolder/*; do
	numb=$((numb+1))	
	perl mapping_pipeline.pl simple_index${numb}.txt test_run $f 10 fragment_map/ test_repeat/135
done

