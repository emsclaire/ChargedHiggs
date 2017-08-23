#include "include/TMVAClassification.h"

TMVAClassification::TMVAClassification(std::string inputFileName, std::string typeBDT) :
  m_inputFileName(inputFileName),
  m_typeBDT(typeBDT) {
    m_outputTrain = "data/train_" + typeBDT +".root";
    m_outputTMVA  = "tmva_" + typeBDT + ".root";
    ReadInputFiles();
    TMVA::Tools::Instance();
}

void TMVAClassification::ReadInputFiles() {

  m_inputFile = new TFile(m_inputFileName.c_str());

  return;

}

void TMVAClassification::SplitTrees() {

  // Checking if the files with the trees exists

  std::ifstream infile(m_outputTrain.c_str());

  if (infile.good()) {

    std::cout << "The file exists. Going to read the trees : " << m_outputTrain.c_str() << std::endl;

    TFile *treeFile = new TFile(m_outputTrain.c_str());
    m_signalTree  =  (TTree*)  treeFile->Get("Signal");
    m_bkgTree     =  (TTree*)  treeFile->Get("Background");
    return;

  }


  TFile *treeFile = new TFile(m_outputTrain.c_str(),"RECREATE");

  TTree* originalTree = (TTree*) m_inputFile->Get("RecoResults");

  std::cout << " --- Creating Signal Tree" << std::endl;

  m_signalTree  = originalTree->CopyTree("CorrectMatch == 1");
  m_signalTree->SetName("Signal");

  std::cout << " --- Creating Background Tree" << std::endl;

  m_bkgTree     = originalTree->CopyTree("CorrectMatch == 0");
  m_bkgTree->SetName("Background");

  std::cout << "Number of events in signal tree: "<< m_signalTree->GetEntries() << std::endl;
  std::cout << "Number of events in background tree: "<< m_bkgTree->GetEntries() << std::endl;

  treeFile->Write();

}


void TMVAClassification::Book() {

  m_outputFile = new TFile(m_outputTMVA.c_str(),"RECREATE");

  m_factory = new TMVA::Factory( "TMVAClassification", m_outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification" );

  // m_factory->AddSpectator("EventNumber");
  // m_factory->AddSpectator("NumberOfJets");
  // m_factory->AddSpectator("NumberBTaggedJets");
  // m_factory->AddSpectator("Index_HcJet");
  // m_factory->AddSpectator("Index_tHcJet");
  // m_factory->AddSpectator("Index_tJet");
  // m_factory->AddSpectator("Index_gJet");
  // m_factory->AddSpectator("CorrectMatch");


  m_factory->AddVariable("Pt_HcJet");
  m_factory->AddVariable("Pt_tHcJet");
  m_factory->AddVariable("Pt_tJet");
  m_factory->AddVariable("Pt_gJet");
  m_factory->AddVariable("Eta_HcJet");
  m_factory->AddVariable("Eta_tHcJet");
  m_factory->AddVariable("Eta_tJet");
  m_factory->AddVariable("Eta_gJet");
  m_factory->AddVariable("Phi_HcJet");
  m_factory->AddVariable("Phi_tHcJet");
  m_factory->AddVariable("Phi_tJet");
  m_factory->AddVariable("Phi_gJet");
  m_factory->AddVariable("E_HcJet");
  m_factory->AddVariable("E_tHcJet");
  m_factory->AddVariable("E_tJet");
  m_factory->AddVariable("E_gJet");
  m_factory->AddVariable("dPtPair_HcJet_tHcJet");
  m_factory->AddVariable("dPtPair_HcJet_tJet");
  m_factory->AddVariable("dPtPair_HcJet_gJet");
  m_factory->AddVariable("dPtPair_tHcJet_tJet");
  m_factory->AddVariable("dPtPair_tHcJet_gJet");
  m_factory->AddVariable("dPtPair_tJet_gJet");
  m_factory->AddVariable("dEtaPair_HcJet_tHcJet");
  m_factory->AddVariable("dEtaPair_HcJet_tJet");
  m_factory->AddVariable("dEtaPair_HcJet_gJet");
  m_factory->AddVariable("dEtaPair_tHcJet_tJet");
  m_factory->AddVariable("dEtaPair_tHcJet_gJet");
  m_factory->AddVariable("dEtaPair_tJet_gJet");
  m_factory->AddVariable("dPhiPair_HcJet_tHcJet");
  m_factory->AddVariable("dPhiPair_HcJet_tJet");
  m_factory->AddVariable("dPhiPair_HcJet_gJet");
  m_factory->AddVariable("dPhiPair_tHcJet_tJet");
  m_factory->AddVariable("dPhiPair_tHcJet_gJet");
  m_factory->AddVariable("dPhiPair_tJet_gJet");
  m_factory->AddVariable("dMPair_HcJet_tHcJet");
  m_factory->AddVariable("dMPair_HcJet_tJet");
  m_factory->AddVariable("dMPair_HcJet_gJet");
  m_factory->AddVariable("dMPair_tHcJet_tJet");
  m_factory->AddVariable("dMPair_tHcJet_gJet");
  m_factory->AddVariable("dMPair_tJet_gJet");
  m_factory->AddVariable("dRPair_HcJet_tHcJet");
  m_factory->AddVariable("dRPair_HcJet_tJet");
  m_factory->AddVariable("dRPair_HcJet_gJet");
  m_factory->AddVariable("dRPair_tHcJet_tJet");
  m_factory->AddVariable("dRPair_tHcJet_gJet");
  m_factory->AddVariable("dRPair_tJet_gJet");
  m_factory->AddVariable("PtPair_HcJet_tHcJet");
  m_factory->AddVariable("PtPair_HcJet_tJet");
  m_factory->AddVariable("PtPair_HcJet_gJet");
  m_factory->AddVariable("PtPair_tHcJet_tJet");
  m_factory->AddVariable("PtPair_tHcJet_gJet");
  m_factory->AddVariable("PtPair_tJet_gJet");
  m_factory->AddVariable("EtaPair_HcJet_tHcJet");
  m_factory->AddVariable("EtaPair_HcJet_tJet");
  m_factory->AddVariable("EtaPair_HcJet_gJet");
  m_factory->AddVariable("EtaPair_tHcJet_tJet");
  m_factory->AddVariable("EtaPair_tHcJet_gJet");
  m_factory->AddVariable("EtaPair_tJet_gJet");
  m_factory->AddVariable("PhiPair_HcJet_tHcJet");
  m_factory->AddVariable("PhiPair_HcJet_tJet");
  m_factory->AddVariable("PhiPair_HcJet_gJet");
  m_factory->AddVariable("PhiPair_tHcJet_tJet");
  m_factory->AddVariable("PhiPair_tHcJet_gJet");
  m_factory->AddVariable("PhiPair_tJet_gJet");
  m_factory->AddVariable("EPair_HcJet_tHcJet");
  m_factory->AddVariable("EPair_HcJet_tJet");
  m_factory->AddVariable("EPair_HcJet_gJet");
  m_factory->AddVariable("EPair_tHcJet_tJet");
  m_factory->AddVariable("EPair_tHcJet_gJet");
  m_factory->AddVariable("EPair_tJet_gJet");
  m_factory->AddVariable("Minv_HcJet_PosLep");
  m_factory->AddVariable("Minv_tHcJet_PosLep");
  m_factory->AddVariable("Minv_tJet_PosLep");
  m_factory->AddVariable("Minv_gJet_PosLep");
  m_factory->AddVariable("Minv_HcJet_NegLep");
  m_factory->AddVariable("Minv_tHcJet_NegLep");
  m_factory->AddVariable("Minv_tJet_NegLep");
  m_factory->AddVariable("Minv_gJet_NegLep");
  m_factory->AddVariable("PtPair_HcJet_PosLep");
  m_factory->AddVariable("PtPair_tHcJet_PosLep");
  m_factory->AddVariable("PtPair_tJet_PosLep");
  m_factory->AddVariable("PtPair_gJet_PosLep");
  m_factory->AddVariable("PtPair_HcJet_NegLep");
  m_factory->AddVariable("PtPair_tHcJet_NegLep");
  m_factory->AddVariable("PtPair_tJet_NegLep");
  m_factory->AddVariable("PtPair_gJet_NegLep");
  m_factory->AddVariable("EtaPair_HcJet_PosLep");
  m_factory->AddVariable("EtaPair_tHcJet_PosLep");
  m_factory->AddVariable("EtaPair_tJet_PosLep");
  m_factory->AddVariable("EtaPair_gJet_PosLep");
  m_factory->AddVariable("EtaPair_HcJet_NegLep");
  m_factory->AddVariable("EtaPair_tHcJet_NegLep");
  m_factory->AddVariable("EtaPair_tJet_NegLep");
  m_factory->AddVariable("EtaPair_gJet_NegLep");
  m_factory->AddVariable("PhiPair_HcJet_PosLep");
  m_factory->AddVariable("PhiPair_tHcJet_PosLep");
  m_factory->AddVariable("PhiPair_tJet_PosLep");
  m_factory->AddVariable("PhiPair_gJet_PosLep");
  m_factory->AddVariable("PhiPair_HcJet_NegLep");
  m_factory->AddVariable("PhiPair_tHcJet_NegLep");
  m_factory->AddVariable("PhiPair_tJet_NegLep");
  m_factory->AddVariable("PhiPair_gJet_NegLep");
  m_factory->AddVariable("EPair_HcJet_PosLep");
  m_factory->AddVariable("EPair_tHcJet_PosLep");
  m_factory->AddVariable("EPair_tJet_PosLep");
  m_factory->AddVariable("EPair_gJet_PosLep");
  m_factory->AddVariable("EPair_HcJet_NegLep");
  m_factory->AddVariable("EPair_tHcJet_NegLep");
  m_factory->AddVariable("EPair_tJet_NegLep");
  m_factory->AddVariable("EPair_gJet_NegLep");


   // global event weights per tree (see below for setting event-wise weights)
   double signalWeight     = 1.0;
   double backgroundWeight = 1.0;
   // You can add an arbitrary number of signal or background trees
   m_factory->AddSignalTree    ( m_signalTree ,     signalWeight );
   m_factory->AddBackgroundTree( m_bkgTree    , backgroundWeight );


  m_factory->BookMethod(  TMVA::Types::kBDT, "BDT_" + m_typeBDT,
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
  SplitTrees();
  Book();
  TrainAndTest();
  return;
}
