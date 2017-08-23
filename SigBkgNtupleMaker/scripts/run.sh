#!/bin/bash

sample=$1

INFILE="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/RecoBDTNtuples/${sample}.root"
OUTFILE="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/SigBkgNtuples/${sample}.root"

bin/NtupleMaker ${INFILE} ${OUTFILE}
