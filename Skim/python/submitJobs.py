import sys
import os

def getListsOfNtuples( file ):
	with open(file) as f:
		return f.read().splitlines()


def makeScript( workdir, ntuplesdir, outdir, filename, is_signal):
	body='''#!/bin/bash

WORKDIR="%(workdir)s"
NTUPLEDIR="%(ntuplesdir)s"
OUTDIR="%(outdir)s"

cd $WORKDIR

source scripts/setup.sh

%(workdir)s/bin/Analyser ${NTUPLEDIR}/%(filename)s ${OUTDIR}/%(filename)s %(signal)s no
'''% {'workdir':workdir,'outdir':outdir,'ntuplesdir':ntuplesdir,'filename':filename,'signal': 'sig' if is_signal else 'bkg'}

	with open(workdir + '/batch/' +filename[:-5]+ '.sh','w+') as f:
		f.write(body)

if __name__ == '__main__':

	sample = sys.argv[1] if len(sys.argv)>1  else 'm300' # Sample to run: m900, m300, ttbb, ttcc, ttlight
	numbers_of_samples_to_run = int(sys.argv[2]) if len(sys.argv)>2 else 5 # Need 35 for full, 250 for full split samples
        is_signal = True if 'm' in sample else False

	ntuple_directory = '/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/Ntuples/pp/' + sample
	work_directory   = '/afs/cern.ch/user/e/eorgill/public/ChargedHiggs/Code/Skim'
	out_directory    = '/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/SkimmedNtuples'

	ntuples = getListsOfNtuples(work_directory + '/data/' + sample + '.list')

	for i in range(numbers_of_samples_to_run):
		filename = ntuples[i]
		makeScript( work_directory, ntuple_directory, out_directory, filename, is_signal )
		submit_cmd = 'bsub -R "pool>3000" -q 8nh -J %s < %s/batch/%s.sh' % (filename[:-5],work_directory,filename[:-5])
		os.chdir(work_directory + '/batch')
		os.system(submit_cmd)
