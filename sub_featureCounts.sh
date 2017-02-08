#!/bin/bash
#Align RNAseq data with genome using featureCounts

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 4
#$ -l virtual_free=1.2G

# ---------------
# Step 1
# Collect inputs
# ---------------

InBam=$(basename $1)
InGff=$(basename $2)
OutDir=$3
Prefix=$4

CurDir=$PWD
WorkDir=$TMPDIR/featureCounts
echo "$WorkDir"
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$1 $InBam
cp $CurDir/$2 $InGff
cp $CurDir/$3 $OutDir


# ---------------
# Step 2
# Run featureCounts
# ---------------

/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/subread-1.5.1-source/bin/featureCounts \
-p -B -M -R -a $InGff -t exon -g "Parent" -o "$Prefix"_featurecounts.txt $InBam

rm $InBam
rm $InGff
mkdir -p $CurDir/$OutDir
cp -r $WorkDir $CurDir/$OutDir/.
