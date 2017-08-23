#!/bin/bash
NTUPLEDIR="/afs/cern.ch/work/e/eorgill/public/2HDM_SignalSamples"
WORKDIR="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/NeutrinoWeighting/Analysis"

cd ${WORKDIR}/batch
while read p; do
  bsub -R "pool>3000" -q 8nh -J job1 < ${WORKDIR}/scripts/runAnalyser.sh $p
done < ${WORKDIR}/data/m900.list
