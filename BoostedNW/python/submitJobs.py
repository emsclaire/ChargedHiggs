import sys
import os

def getListsOfNtuples( file ):
	with open(file) as f:
		return f.read().splitlines()


def makeScript( workdir, ntuplesdir, outdir, filename, newfile, sample):
	body='''#!/bin/bash

WORKDIR="%(workdir)s"
NTUPLEDIR="%(ntuplesdir)s"
OUTDIR="%(outdir)s"

cd $WORKDIR

source scripts/setup.sh

%(workdir)s/bin/BNWReco ${NTUPLEDIR}/%(filename)s ${OUTDIR}/%(newfilename)s %(sample)s
'''% {'workdir':workdir,'outdir':outdir,'ntuplesdir':ntuplesdir,'filename':filename,'newfilename':newfile,'sample': sample}

	with open(workdir + '/batch/' +newfilename[:-5]+ '.sh','w+') as f:
		f.write(body)

if __name__ == '__main__':

	sample = sys.argv[1] if len(sys.argv)>1  else 'm300' # Sample to run: m900 or m300
	numbers_of_samples_to_run = int(sys.argv[2]) if len(sys.argv)>2 else 10 # Need 35 for full samples
        is_signal = True if "m" in sample else False
        typeBDT = sys.argv[3] if not "m" in sample else sample

	ntuple_directory = '/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/SkimmedNtuples'
	work_directory   = '/afs/cern.ch/user/e/eorgill/public/ChargedHiggs/Code/BoostedNW'
	out_directory    = '/afs/cern.ch/work/e/eorgill/public/ChargedHiggs/RecoBDTNtuples'

	ntuples = getListsOfNtuples(work_directory + '/data/' + sample + '.list')

	for i in range(numbers_of_samples_to_run):
		filename = ntuples[i]
                newfilename = filename if sample == typeBDT else typeBDT +"_"+filename
		makeScript( work_directory, ntuple_directory, out_directory, filename, newfilename, typeBDT )
		submit_cmd = 'bsub -R "pool>3000" -q 8nh -J %s < %s/batch/%s.sh' % (filename[:-5],work_directory,newfilename[:-5])
		os.chdir(work_directory + '/batch')
		os.system(submit_cmd)
