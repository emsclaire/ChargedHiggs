#!/bin/bash

sample=$1
typeBDT=$2
bdtcut=$3

INFILE="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/RecoBDTNtuples/${sample}.root"
OUTFILE="/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/TRExFitter/InputHists/geq4jgeq4b/${sample}_BDTCUT_${bdtcut}.root"

bin/MakeHists ${INFILE} ${OUTFILE} ${typeBDT} ${bdtcut}
