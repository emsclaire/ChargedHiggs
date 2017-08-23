#!/bin/bash

WORKDIR=/afs/cern.ch/user/e/eorgill/public/ChargedHiggs/Code/Skim
INDIR=/afs/cern.ch/user/e/eorgill/work/public/ChargedHiggs/Ntuples
OUTDIR=/afs/cern.ch//work/e/eorgill/public/ChargedHiggs/SplittedNtuples

while read file ; do
echo "procesing file: $file "


outFilename=`echo $file | cut -d "." -f 1`



echo "splitting the samples"
$WORKDIR/bin/treeSplit $INDIR/$file 1000 $OUTDIR/$outFilename


done < $WORKDIR/data/allSignal.list
