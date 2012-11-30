#! /bin/bash 


# just a small shell script, that systematically evaluates all
# datasets and then providse the 4 plots.

# We have either ML or parsimony and either hard or easy datasets.



for ml in  True  # False 
do 
    for hardness in easy hard 
    do 
	outputFile=$hardness.evaluated
	rm -f  $outputFile
	for someFile in $(ls ../data/simulated/$hardness/info.*)
	do 
	    echo $someFile
	    theid=$(basename $someFile .txt | cut -f 2,3,4,5 -d '.') 
	    dataset=$(basename $someFile .txt | cut -f 2 -d '.')	
	    ./evaluateADataset.py ../data/origData/$dataset.tre  ../data/origData/$dataset.aln  ../data/simulated/easy/mislabeledTaxonomy.$theid.tax $ml ../data/simulated/easy/info.$theid.txt >> $outputFile
	done

	./plotROC.R $hardness.$ml $outputFile > /dev/null 
    done 
done 
