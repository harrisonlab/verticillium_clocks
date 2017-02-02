#!/bin/bash
#Align RNAseq data with genome using STAR

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 16
#$ -l virtuak_free=1.2G
#$ -l

# ---------------
# Step 1
# Collect inputs
# ---------------

InGenome=$(basname $1)
InGff=$(basname $2)
InReadF=$(basname $3)
InReadR=$(basname $4)
GenomeDir=$(basname $5)
OutDir=$6

CurDir=$PWD
WorkDir=$TMPDIR/star
mkdir -p $WorkDir
cd $WorkDir

cp $CurDir/$1 $InGenome
cp $CurDir/$2 $InGff
cp $CurDir/$3 $InReadF
cp $CurDir/$4 $InReadR
cp $CurDir/$5 $GenomeDir
cp $CurDir/$6 $OutDir

# ---------------
# Step 2
# Create the Index File
# ---------------

ParentFeature="Parent"
# ParentFeature="ID"

./STAR \
--runMode genomeGenerate \
--genomeDir $GenomeDir \
--genomeFastaFiles $InGenome \
--sjdbGTFtagExonParentTranscript $ParentFeature \
--sjdbGTFfile $InGff
--runThreadN 16 \
--sjdbOverhang ReadLength-1

# ---------------
# Step 2=3
# Run STAR
# ---------------

./STAR \
--genomeDir $GenomeDir \
--outFileNamePrefix $OutDir \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype SAM \
--runThreadN 16


rm -r $GenomeDir
rm $InGenome
rm $InGff
rm $InReadF
rm $InReadR
mkdir -p $CurDir/$OutDir
cp -r $WorkDir $CurDir/$OutDir/.
