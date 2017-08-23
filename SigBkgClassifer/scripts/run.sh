#!/bin/bash

sample=$1
region=$2

SIGFILE="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/SigBkgNtuples/${sample}.root"
BKGFILE="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/SigBkgNtuples/${sample}_background.root"

bin/BDTTrainer ${SIGFILE} ${BKGFILE} ${sample} ${region}
