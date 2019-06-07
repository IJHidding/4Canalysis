#!/bin/bash
#SBATCH --job-name=DataAnalysis
#SBATCH --time=11:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=30gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L



#input should be all the files you want analyzed. folder containing input files

fastqfiles=/groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/fastqdir
primerseqfile=/groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/Primersequencefile


module load cutadapt
module load Bowtie2
#echo "demultiplexing viewpoint"

#demultiplexing on all primers included in the primerseqfile, up to 8.
for f in /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/${fastqfiles}/*.fq; do
	#echo $(awk 'NR==1 {print $2}' ${primersequencefile})
	cutadapt -O 13 -g $(awk 'NR==1 {print $2}' ${primerseqfile}) -g $(awk 'NR==2 {print $2}' ${primerseqfile}) -g $(awk 'NR==3 {print $2}' ${primerseqfile}) -g $(awk 'NR==4 {print $2}' ${primerseqfile}) -g $(awk 'NR==5 {print $2}' ${primerseqfile}) -g $(awk 'NR==6 {print $2}' ${primerseqfile}) -g $(awk 'NR==7 {print $2}' ${primerseqfile}) -g $(awk 'NR==8 {print $2}' ${primerseqfile}) -o '{name}' $f	
	filename=$(basename -- "$f")
	InputnameNoExtention="${filename%.*}"

# making some tmp storage folders
mkdir ./tmp
mkdir ./unalign
mkdir ./align
mkdir ./finished

#renaming cutadapt output
	echo "Renaming files"
	mv 1 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==1 {print $1}' ${primerseqfile}).fq
	mv 2 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==2 {print $1}' ${primerseqfile}).fq
	mv 3 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==3 {print $1}' ${primerseqfile}).fq
	mv 4 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==4 {print $1}' ${primerseqfile}).fq
	mv 5 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==5 {print $1}' ${primerseqfile}).fq
	mv 6 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==6 {print $1}' ${primerseqfile}).fq
	mv 7 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==7 {print $1}' ${primerseqfile}).fq
	mv 8 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/${InputnameNoExtention}_$(awk 'NR==8 {print $1}' ${primerseqfile}).fq
	
done
	

#Alignment using bowtie2 with a preindexed Human_v37 genome stores the unaligned reads for later.
echo "First alignment"
for file in /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/*; do
	basefile=$(basename -- "$file")
	echo $basefile   
	bowtie2 -x ../4C/data/bowtie2index/Human_v37 $file -S /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/finished/${basefile}.sam --un /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/unalign/unaligned_${basefile}
done
	
# split unaligned data on restriction sites, added after the primers in the sequence file and renaming each file.
echo "splitting unaligned"
for unaligned in /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/unalign/*.fq; do
	basealigned=$(basename -- "$unaligned")
	cutadapt -g $(awk 'NR==10 {print $2}' ${primerseqfile}) -g $(awk 'NR==9 {print $2}' ${primerseqfile}) -o '{name}' $unaligned
	mv 1 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/align/$(awk 'NR==9 {print $2}' ${primerseqfile})_rightside_${basealigned}
	mv 2 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/align/$(awk 'NR==10 {print $2}' ${primerseqfile})_rightside_${basealigned}
	cutadapt -a $(awk 'NR==10 {print $2}' ${primerseqfile}) -a $(awk 'NR==9 {print $2}' ${primerseqfile}) -o '{name}' $unaligned
	mv 1 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/align/$(awk 'NR==9 {print $2}' ${primerseqfile})_leftside_${basealigned}
        mv 2 /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/align/$(awk 'NR==10 {print $2}' ${primerseqfile})_leftside_${basealigned}
done

#Attemping to align the split files again.
echo "Second alignment"
for second in /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/align/*; do
	basesecond=$(basename -- "$second")
	echo $basesecond
	bowtie2 -x ../4C/data/bowtie2index/Human_v37 $second -S /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/finished/${basesecond}.sam --un /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/unalign2/unaligned_${basesecond}
done

module load SAMtools

#Move samfiles into Bamfiles.
for f in /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/finished/*.sam; do
        filename=$(basename -- "$f")
        InputnameNoExtention="${filename%.*}"
        samtools view -b $f > /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/bamstorage/${InputnameNoExtention}.bam
done


#Merge matching Samfiles to Bamfiles based on the viewpoint FQ files.
for f in /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/tmp/*.fq; do
        filename=$(basename -- "$f")
        InputnameNoExtention="${filename%.*}"
        OGfile=$InputnameNoExtention
	echo $OGfile
	samtools merge /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/mergedbam/${OGfile}.bam /groups/umcg-gdio/tmp04/umcg-ihidding/181212_M01785_0241_000000000-C6Y44_1812_4C_CM/bamstorage/*${OGfile}*
done


