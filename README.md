******************************************************************************************************************

Charged Higgs Pheno Study

******************************************************************************************************************

Name: Emily Orgill

Email: emily.orgill@postgrad.manchester.ac.uk

Date: 21 July 2017

******************************************************************************************************************

Description:

To place limits on ttHc -> bbbbllvv cross-section

Using MadGraph5+Pythia6+Delphes samples (875000 events per sample) (updating to Pythia8)

Signal: Charged Higgs samples produced with 2HDMC model at mass = {300,400,500,600,700,800,900} GeV

Background: ttbb, ttcc and tt+l

******************************************************************************************************************

Code

******************************************************************************************************************

Skim:

(obtain TLorentzVectors, match jets to b-partons and apply cuts etc.)

    python python/submitJobs.sh [sample] [njobs]

    Input:  data/[sample].list listing Delphes root file names
            [sample] = m300, m400, m500, m600, m700, m800, m900, ttbb, ttcc, ttlight
            [njobs] = 35
    Output: SkimmedNtuples

******************************************************************************************************************

RecoTrainingNtupleMaker:

(create input to RecoBDT training)

    ./scripts/run.sh [inputfilelist] [outputfile]

    Input:  [inputfilelist] = [sample]skim.txt (in data dir)
            [sample] = m300, m400, m500, m600, m700, m800, m900
    Output: PermutationNtuples
            [outputfile] = [sample].root

******************************************************************************************************************

BoostedNW -> BDTTrainer:

(train reconstruction BDT for choosing jet permutation with correct source)

    ./bin/BDTTrainer /PATH/TO/PermutationNtuples/[sample].root [sample]

    Input:  PermutationNtuples/[sample].root
            [sample] = m300, m400, m500, m600, m700, m800, m900
    Output: tmva_[sample].root
            weights/TMVAClassification_BDT_[sample].weights.xml

******************************************************************************************************************
Boosted NW -> BNWReco:

(contains NW solution for permutation chosen by BDT)

    python python/submitJobs.sh [sample] [njobs] [bdt]

    Input:  SkimmedNtuples, weights/TMVAClassification_BDT_[sample].weights.xml
            [bdt] = m300, m400, m500, m600, m700, m800, m900
            [sample] = m300, m400, m500, m600, m700, m800, m900, ttbb, ttcc, ttlight
    Output: RecoBDTNtuples

******************************************************************************************************************

SigBkgNtupleMaker:

(create ntuples to feed to classification BDT)

    cd /PATH/TO/RecoBDTNtuples

    hadd m[mass].root run-*_m-[mass]_tb-2_ebeam-6500_delphes_events.root
    hadd m[mass]_ttbb.root m[mass]_ttbb_run-*_delphes_events.root
    hadd m[mass]_ttcc.root m[mass]_ttcc_run-*_delphes_events.root
    hadd m[mass]_ttlight.root m[mass]_ttlight_run-*_delphes_events.root

    hadd m[mass]_background.root m[mass]_ttbb.root m[mass]_ttcc.root m[mass]_ttlight.root
       [mass] = 300, 400, 500, 600, 700, 800, 900

    cd /PATH/TO/SigBkgNtupleMaker

    ./scripts/run.sh [sample]  (if still requires a BDT cut argument, remove this)

    Input:  RecoBDTNtuples
            [sample] = m[mass], m[mass]_background
    Output: SigBkgNtuples

******************************************************************************************************************

SigBkgClassifier -> BDTTrainer

(train signal vs background classification BDT)

    ./scripts/run.sh [sample] [region]

    Input:  SigBkgNtuples
            [sample] = m300, m400, m500, m600, m700, m800, m900
            [region] = inclusive
    Output: data/train_[sample].root
            weights/TMVAClassification_BDT_[sample]_[region].weights.xml

******************************************************************************************************************

SigBkgClassifier -> LimitsHistMaker:

(make histograms for input to fit/exclusion)

Edit runhistmaker.sh to save to either geq4jgeq4b or geq4jeq2b

And edit LimitsHistMaker.cxx (line X) to n_btags >= 4 or n_btags == 2 respectively (plan to automate this soon)

    ./scripts/runhistmaker.sh [sample] [typeBDT] [BDTCut]

    Input:  RecoBDTNtuples/[sample].root
            weights/TMVAClassification_BDT_[typeBDT].weights.xml
            [sample] = [mass], [mass]_ttbb, [mass]_ttcc, [mass]_ttlight
            [typeBDT] = [mass]_inclusive, [mass]_inclusive
            [mass] = m300, m400, m500, m600, m700, m800, m900
    Output: InputHists/geq4jgeq4b or InputHists/geq4jeq2b
