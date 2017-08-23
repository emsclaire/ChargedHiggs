#!/bin/bash

inputfile=$1
outputfile=$2

LISTFILE="/afs/cern.ch/user/e/eorgill/public/ChargedHiggs/Code/RecoTrainingNtupleMaker/data/${inputfile}"
OUTFILE="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/PermutationNtuples/${outputfile}"

bin/NtupleMaker "${LISTFILE}" "${OUTFILE}"
