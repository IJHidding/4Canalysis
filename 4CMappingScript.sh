#!/bin/bash

##Lists the current location to find all files, can be used in an install function
####Preset this to ensure proper working
currentlocation=/Path/to/install/location
Jaspardatafile=/Path/to/JasparDatafile
JasparMatrixInfo=/Path/to/Jasparmatrix/
Rlibrary=/Path/To/R/Library

#Help function, gives information on the used parameters and usability of the program
###################################HELP AND INFORMATION########################################################

function showHelp() {
        #
        # Display commandline help on STDOUT.
        #
	cat <<EOH
===============================================================================================================
Pipeline to demultiplex, align and analyze 4C data.
Usage:
	sh $(basename $0) OPTIONS
Options:
	-h	Show this help.
	Required:
	-i	inputfile (See readme for details)
	Optional:
	-w	Specify the size of the running window (default=300)
	-l	Specify the steps taken in the running window (default=50)
	-r	Runs the R installer to install the peakC package and required libraries
===============================================================================================================
EOH
	trap - EXIT
        exit 0
}
################################################################################################################
#Checks if the folders exists and creates them if they dont
#################################################################################################################
function Prepwork() {
	if [ ! -d "${currentlocation}/4C_Output" ];
	then
		mkdir ${currentlocation}/4C_Output
	fi

	output=${currentlocation}/4C_Output
	if [ ! -d "${output}/AnalysisOutput" ];
	then
	        mkdir ${output}/AnalysisOutput
	fi

	Jasparoutput=${output}/AnalysisOutput

	PathToReferenceGenome=${currentlocation}/

	data=${currentlocation}/inputfiles

	
	#Determines the length of the reads in the file.
	reads=$(awk 'NR==2' ${data}/*  | wc -m)
	lengthofread=$((${reads} -1))

	if [ ! -d "${output}/test_run/splitted_seq/" ]; 
	then 
		mkdir ${output}/test_run
		mkdir ${output}/test_run/splitted_seq
	fi

	#Handling the input file into a cutadapt file and a mapping pipeline file
	## finding the used restriction sites to generate a fragment map.
	RES1=$(cut -f4 $input | tail -1)
	RES2=$(cut -f5 $input | tail -1) 
}
##########################################################################################################################
#The cutadapt function
#Runs the cutadapt program based on the amount of primers included in the input files. 
##########################################################################################################################
function Demultiplexer() {
	number=0
	#Demultiplexing fastq files using cutadapt 
	#Names the files after the primer names.
	for inputfile in  ${data}/*.f*; do
		number=$((${number} +1))
		for primerseq in $(cut -f2 $input); do
			primerseq=$(echo $primerseq | sed 's/GATC//g')
			Primerlength=$(echo $primerseq| wc -m) 
			cutadapt -O 13 -g $primerseq -o '{name}' $inputfile 
			PrimerName=$(grep ${primerseq} $input | cut -f1)
			mv 1 ${output}/test_run/splitted_seq/${PrimerName}${number}.fastq
			cp ${output}/test_run/splitted_seq/${PrimerName}${number}.fastq ${output}/
		done
	done
	echo "Demultiplex check, so far so good."
}
##########################################################################################################################
#Aligning and alignment prep.
##########################################################################################################################
function AllTheOtherStuff() {
	#Calculates readlength without primers but including the restriction site.
	PrimerLength=$(awk 'NR==1' $input | cut -f2 | wc -m)
	ActualPrimerLength=$(echo $((${lengthofread}-${Primerlength} + 4)))
	  
	#prep fragment map
	#check if fragment map exists
	if [ ! -d "${currentlocation}/fragment_map" ];
	then
		perl ${currentlocation}/generate_fragment_map.pl ${PathToReferenceGenome}ucsc.hg19.fa ${RES1} ${RES2} ${currentlocation}/fragment_map/
		echo "Fragment map has been succesfully made!"
	else
		echo "Fragment map already exists!"
	fi


	#check if repeats are already made
	if [ ! -d "${currentlocation}/test_repeat/${ActualPrimerLength}" ];
	then 
		perl ${currentlocation}/getRepeats.pl ${currentlocation}/fragment_map/ ${RES1} ${ActualPrimerLength} ${PathToReferenceGenome}ucsc.hg19.fa ${currentlocation}/test_repeat/
		echo "Repeats for the ${ActualPrimerLength} are made!"
	else
		echo "Repeats for the ${ActualPrimerLength} length already exist!"
	fi

	#Creates a file for each fastq file that numbers the input files for analysis.
	NumberofInputfiles=$(ls ${currentlocation}/inputfiles/ | wc -l) 
	for Amount in $(seq ${NumberofInputfiles}); do
		awk -v var="$Amount" '$1=$1var' ${input} > ${input}${Amount}.txt
		sed -i 's/ /\t/g' ${input}${Amount}.txt
	done
	echo "Start of alignment!"

	#Runs the mapping pipeline to map the fastq data.
	numberoffiles=0
	for f in ${currentlocation}/inputfiles/*; do
		numberoffiles=$((numberoffiles+1))
		perl ${currentlocation}/mapping_pipeline.pl ${input}${numberoffiles}.txt ${output}/test_run $f 10 ${currentlocation}/fragment_map/ ${currentlocation}/test_repeat/${ActualPrimerLength}	
		rm ${input}${numberoffiles}.txt
	done
	rm ${currentlocation}/unknown

}



##################################################################################################################
#Function that runs the Rscript and does the transcription factor prioritization.
#Does the peakC analysis and performs the prioritization per viewpoint.
#
##################################################################################################################
function RandAnalysis() {
	echo "Rscriptstart"
	#PeakC
	###Initializing R
	export R_LIBS=${Rlibrary}
	for viewpoint in $(cut -f1 ${input}); do

		Pattern=${viewpoint}
		viewpointnumber=$(grep ${viewpoint} ${input} | cut -f7)
		viewpointfile=${output}/test_run/filtered_file/*${viewpoint}*
		#creates an output file for each viewpoint
		outputfile=${output}/${viewpoint}
		if [ ! -d "${outputfile}" ];
		then
			mkdir ${outputfile}
		fi
		Rscript ${currentlocation}/PeakC.R $viewpoint $viewpointnumber $outputfile/ ${output}/test_run/filtered_file/ ${Pattern} ${outputfile}/${Pattern}.txt
		mv $outputfile/Plot.jpg $outputfile/${Pattern}.jpg
	

		Jasparinput=${outputfile}/${viewpoint}.txt
		NameofInputFile=${inputname}
		filename=$(basename -- "$Jasparinput")
		InputnameNoExtention="${filename%.*}"
		chromosomenumber=$(grep $InputnameNoExtention ${input} | cut -f6)
		NumberofPeaks=0
		sed -i '/x/d' $Jasparinput
		if [ -s $Jasparinput ]; 
		then 
			coordinate2=0
			for coordinate in $(cat $Jasparinput); do
	 			coordinate1=$coordinate
				ComparedCoordinate=$(($coordinate1 - $coordinate2))
				if [ "${ComparedCoordinate}" -ge "30000" ];
				then
	                		NumberofPeaks=$(($NumberofPeaks+1))
	        		fi
	        		echo $coordinate >> ${Jasparoutput}/${InputnameNoExtention}${NumberofPeaks}.txt
	        		coordinate2=$coordinate1
			done
	
			#for numbers in peak file take highest and lowest
			#Create a running window of 300 bp

			for peak in ${Jasparoutput}/${InputnameNoExtention}*; do
				echo $peak
		        	peakername=$(basename -- "$peak")
		        	PeakName="${peakername%.*}"
		        	start=$(cat $peak | head -1)
		        	end=$(cat $peak | tail -1)
		        	number1=$start
		        	while [ "$number1" -lt "$end" ]; do
		        	        number2=$(($number1+${Window}))
		        	        echo "${chromosomenumber}       ${number1}      ${number2}      ${PeakName}" >> ${Jasparoutput}/${InputnameNoExtention}.bed
		        	        number1=$((${number1}+${LengthofWindow}))
		        	done
			done

			#Compare Jaspar database using bedtools
			for f in ${Jaspardatafile}MA*; do
			        MatrixName=$(basename -- "$f")
			        MatrixNameNoExtention="${MatrixName%.*}"
				sed -i 's/ \+ /\t/g' ${Jasparoutput}/${InputnameNoExtention}.bed
			        bedtools intersect -wa -a ${Jasparoutput}/${InputnameNoExtention}.bed -wb -b $f > ${Jasparoutput}/${MatrixNameNoExtention}.tmp
			        TranscriptionFactor=$(awk 'NR==5' ${JasparMatrixInfo}${MatrixNameNoExtention}.jaspar)
			        sed -i "s/$/\t${TranscriptionFactor}/" ${Jasparoutput}/${MatrixNameNoExtention}.tmp
			done
			cat ${Jasparoutput}/*.tmp >> ${Jasparoutput}/${InputnameNoExtention}_output.bed
			sed -i '/hg38/d' ${Jasparoutput}/${InputnameNoExtention}_output.bed
			sort ${Jasparoutput}/${InputnameNoExtention}_output.bed > tmp && mv tmp ${Jasparoutput}/${InputnameNoExtention}_output.bed
			rm ${Jasparoutput}/*.tmp
			rm ${Jasparoutput}/*.txt
		fi	
	done
}
###########################################################################################################
#Installs Rlibrary, and packages
###########################################################################################################
function Rprep() {
        Rscript -e "install.packages('devtools', '${Rlibrary}', 'http://ftp.ussg.iu.edu/CRAN')"
        Rscript -e "install.packages('isotone', '${Rlibrary}', 'http://ftp.ussg.iu.edu/CRAN')"
        Rscript -e "install.packages('caTools', '${Rlibrary}', 'http://ftp.ussg.iu.edu/CRAN')"
        Rscript -e "devtools::install_github('deWitLab/peakC')"
}
###########################################################################################################



#main
while getopts "i:h:wlr" opt;
do
        case $opt in h)showHelp;; i)input="${OPTARG}";; w)Window="${OPTARG}";; l)LengthofWindow="${OPTARG}";; r)Rprep;;
        esac
done
if [[ -z "${input:-}" ]]; then showHelp ; echo "No input is given" ; fi
if [[ -z "${Window:-}" ]]; then Window="300" ; fi
if [[ -z "${LengthofWindow:-}" ]]; then LengthofWindow="50" ; fi





#Running the program. 

Prepwork
Demultiplexer
AllTheOtherStuff
RandAnalysis


echo "All done, have a nice day!"

