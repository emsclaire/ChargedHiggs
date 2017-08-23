#!/bin/bash
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/afs/cern.ch/work/r/rnaranjo/public/ChargedHiggs/Delphes-3.3.2/"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/afs/cern.ch/work/r/rnaranjo/public/ChargedHiggs/ExRootAnalysis/"
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase
source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh
lsetup root