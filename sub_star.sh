#!/bin/bash
#Align RNAseq data with genome using STAR

#$ -S /bin/bash
#$ -cwd
#$ -pe smp 8
#$ -l virtual_free=1.2G

# ---------------
# Step 1
# Collect inputs
# ---------------
echo "Hello"
InGenome=$(basename $1)
InGff=$(basename $2)
InReadF=$(basename $3)
InReadR=$(basename $4)
# genomeDir=$(basname $5)
OutDir=$5

CurDir=$PWD
WorkDir=$TMPDIR/star
genomeDir=$WorkDir/index
echo "$WorkDir"
echo "$genomeDir"
mkdir -p $genomeDir
ls $genomeDir
cd $WorkDir

cp $CurDir/$1 $InGenome
cp $CurDir/$2 $InGff
cp $CurDir/$3 $InReadF
cp $CurDir/$4 $InReadR
# cp $CurDir/$5 $GenomeDir
cp $CurDir/$5 $OutDir



# ---------------
# Step 2
# Create the Index File
# ---------------
echo "Building index file"
ParentFeature="Parent"
# ParentFeature="ID"

/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/STAR \
--runMode genomeGenerate \
--genomeDir $genomeDir \
--genomeFastaFiles $InGenome \
--sjdbGTFtagExonParentTranscript $ParentFeature \
--sjdbGTFfile $InGff \
--runThreadN 8 \
--sjdbOverhang 99

# ---------------
# Step 2=3
# Run STAR
# ---------------

echo "Aligning RNAseq reads"

/home/groups/harrisonlab/project_files/verticillium_dahliae/clocks/STAR \
--genomeDir $genomeDir \
--outFileNamePrefix star_aligment \
--readFilesCommand zcat \
--readFilesIn $InReadF $InReadR \
--outSAMtype SAM \
--runThreadN 8


rm -r $genomeDir
rm $InGenome
rm $InGff
rm $InReadF
rm $InReadR
mkdir -p $CurDir/$OutDir
cp -r $WorkDir $CurDir/$OutDir/.
