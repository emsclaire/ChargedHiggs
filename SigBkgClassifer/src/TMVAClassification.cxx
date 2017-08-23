#include "include/TMVAClassification.h"

TMVAClassification::TMVAClassification(std::string inputSigFileName, std::string inputBkgFileName, std::string typeBDT, std::string region) :
  m_inputSigFileName(inputSigFileName), m_inputBkgFileName(inputBkgFileName),
  m_typeBDT(typeBDT), m_region(region) {
    m_outputTrain = "data/train_" + typeBDT +".root";
    m_outputTMVA  = "tmva_" + typeBDT + "_"+region+".root";
    ReadInputFiles();
    TMVA::Tools::Instance();
}

void TMVAClassification::ReadInputFiles() {
  m_inputSigFile = new TFile(m_inputSigFileName.c_str());
  m_inputBkgFile = new TFile(m_inputBkgFileName.c_str());
  return;
}

void TMVAClassification::GetSigBkgTrees() {

  TFile *treeFile = new TFile(m_outputTrain.c_str(),"RECREATE");

  m_sigTree = (TTree*)m_inputSigFile->Get("RecoResults");
  m_sigTree->SetName("Signal");
  m_bkgTree = (TTree*)m_inputBkgFile->Get("RecoResults");
  m_bkgTree->SetName("Background");

  std::cout << "Number of events in signal tree: "<< m_sigTree->GetEntries() << std::endl;
  std::cout << "Number of events in background tree: "<< m_bkgTree->GetEntries() << std::endl;

  treeFile->Write();

}


void TMVAClassification::Book() {

  m_outputFile = new TFile(m_outputTMVA.c_str(),"RECREATE");

  m_factory = new TMVA::Factory( "TMVAClassification", m_outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );

  // m_factory->AddSpectator("EventNumber");
  // m_factory->AddVariable("NumberOfJets");
  m_factory->AddVariable("MaxBDTWeight");
  // m_factory->AddVariable("Minv_bLeading_ClosestTop");
  m_factory->AddVariable("Minv_bH_ClosestTop");
  m_factory->AddVariable("Minv_minDRJets");
  // m_factory->AddVariable("Ht");
  m_factory->AddVariable("MinMinv_poslep_b");
  m_factory->AddVariable("MinMinv_neglep_b");
  // m_factory->AddVariable("Pt_LeadingJet");
  m_factory->AddVariable("Pt_bH");
  m_factory->AddVariable("Minv_Higgs");
  m_factory->AddVariable("DE_NonTopJets");
  m_factory->AddVariable("DR_Tops");
  m_factory->AddVariable("DR_bH_tH");
  m_factory->AddVariable("DR_bLeading_TopLeading");
  m_factory->AddVariable("Centrality");
  m_factory->AddVariable("DR_tH_lt");
  m_factory->AddVariable("DR_lep_maxDR_t_bH");
  m_factory->AddVariable("costheta_lepton_jet");
  m_factory->AddVariable("EtabH");
  m_factory->AddVariable("EtabtH");
  m_factory->AddVariable("DeltaEta_bg_bt");
  m_factory->AddVariable("DeltaPhi_btH_bH");
  m_factory->AddVariable("DeltaPhi_bt_bH");
  m_factory->AddVariable("DR_bH_btH");
  m_factory->AddVariable("DR_bg_bt");
  m_factory->AddVariable("DeltaPhi_bt_btH");
  m_factory->AddVariable("DR_bg_bH");
  m_factory->AddVariable("DeltaPhi_poslep_neglep");
  m_factory->AddVariable("DeltaPhi_ltH_bH");
  m_factory->AddVariable("DR_ltH_bH");
  m_factory->AddVariable("DeltaPhi_ltH_bt");
  m_factory->AddVariable("DR_ltH_bt");
  // global event weights per tree (see below for setting event-wise weights)
  double signalWeight     = 1.0;
  double backgroundWeight = 1.0;
  // You can add an arbitrary number of signal or background trees
  m_factory->AddSignalTree(m_sigTree, signalWeight);
  m_factory->AddBackgroundTree(m_bkgTree, backgroundWeight);
  // Apply additional cuts on the signal and background samples (can be different)

  TCut mycuts = "NumberOfJets >= 4 && NumberBTaggedJets >=2"; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  //if (m_region.compare("inclusive"))        mycuts = "NumberOfJets >= 4 && NumberBTaggedJets >=2";
  //if (m_region.compare("eq4jet_eq2bjet"))   mycuts = "NumberOfJets == 4 && NumberBTaggedJets ==2";
  //if (m_region.compare("eq4jet_leq2bjet"))  mycuts = "NumberOfJets == 4 && NumberBTaggedJets >=2";
  //if (m_region.compare("leq4jet_leq3bjet")) mycuts = "NumberOfJets >= 4 && NumberBTaggedJets >=3";
  //if (m_region.compare("eq4jet_eq3bjet"))   mycuts = "NumberOfJets == 4 && NumberBTaggedJets ==3";
  //if (m_region.compare("leq4jet_leq4bjet")) mycuts = "NumberOfJets >= 4 && NumberBTaggedJets >=4";

  // Tell the factory how to use the training and testing events
  //
  // If no numbers of events are given, half of the events in the tree are used
  // for training, and the other half for testing:
  //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  // To also specify the number of testing events, use:
  //    factory->PrepareTrainingAndTestTree( mycut,
  //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
  m_factory->PrepareTrainingAndTestTree( mycuts, mycuts,
                                          "SplitMode=Random:NormMode=NumEvents:!V" );
  m_factory->BookMethod(TMVA::Types::kBDT, "BDT_" + m_typeBDT + "_" + m_region,
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

}


void TMVAClassification::TrainAndTest() {
  // Train MVAs using the set of training events
   m_factory->TrainAllMethods();
   // Evaluate all MVAs using the set of test events
   m_factory->TestAllMethods();
   // Evaluate and compare performance of all configured MVAs
   m_factory->EvaluateAllMethods();
   // --------------------------------------------------------------
   // Save the output
   m_outputFile->Close();
   std::cout << "==> Wrote root file: " << m_outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete m_factory;
}

void TMVAClassification::Run() {
  GetSigBkgTrees();
  Book();
  TrainAndTest();
  return;
}
