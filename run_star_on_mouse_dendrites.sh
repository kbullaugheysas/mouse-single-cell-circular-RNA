#!/bin/bash

sample=$1
starVersion="STAR-2.7.2b"
sharedMem="NoSharedMemory"
STAR="/home/bullaugh/tools/bin/$starVersion"
baseOut="realignment/mouse-dendrites"

if [ -z "$sample" ] ; then
  echo "Must specify sample" >&2
  exit 1
fi

if [ ! -x $STAR ] ; then
  echo "Cannot find STAR: ${STAR}" >&2
  exit 1
fi

genomeDir="/kimdata/bullaugh/resources/private-indexes/star272.mm10.75.fusion"
trimdir="/home/bullaugh/analysis/circular-RNA/mouse-dendrites/filtered-fastq"
outdir="$baseOut/${sample}"
if [ ! -d $trimdir ] ; then
  echo "Expecting ${trimdir} directory to exist" >&2
  exit 1
fi
if [ ! -d $outdir ] ; then
  mkdir -p $outdir
fi
inFiles="${trimdir}/${sample}_1.fq.gz ${trimdir}/${sample}_2.fq.gz"

for i in $inFiles; do
  if [ ! -f $i ] ; then
    echo "missing input fastq file $i" >&2
    exit 1
  fi
done
echo "found input files: ${inFiles}"

cmd="$STAR \
    --outFilterScoreMinOverLread 0 \
    --outFilterMatchNminOverLread 0.4 \
    --twopassMode Basic \
    --chimSegmentMin 10 \
    --chimJunctionOverhangMin 12 \
    --alignSJDBoverhangMin 10 \
    --chimSegmentReadGapMax 3 \
    --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimMultimapScoreRange 10 \
    --chimMultimapNmax 10 \
    --chimNonchimScoreDropMin 20 \
    --chimOutJunctionFormat 1 \
    --runThreadN 8 \
    --genomeDir ${genomeDir} \
    --readFilesCommand zcat \
    --readFilesIn $inFiles \
    --genomeLoad $sharedMem \
    --outSAMtype BAM Unsorted \
    --outFileNamePrefix ${outdir}/"

echo $cmd
eval $cmd

# END
