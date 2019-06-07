Jaspardatafile=/groups/umcg-gdio/tmp04/umcg-ihidding/IWAN/Jasparbed/

for f in ${Jaspardatafile}MA*; do
	echo $f >> TFSoutputinrgion.txt
	bedtools intersect -wb -b $1 -wa -a $f >> TFSoutputinrgion.txt
done
