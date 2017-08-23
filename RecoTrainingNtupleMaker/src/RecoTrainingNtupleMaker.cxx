#include "include/RecoTrainingNtupleMaker.h"

RecoTrainingNtupleMaker::RecoTrainingNtupleMaker(std::string inputListFileName, std::string outputFileName) : m_inputListFileName(inputListFileName), m_outputFileName(outputFileName) {
  ReadInputFiles();
}

void RecoTrainingNtupleMaker::ReadInputFiles() {

  std::string line;
  std::ifstream inputFile;
  inputFile.open(m_inputListFileName);

  while(std::getline(inputFile,line)) {
    if (line.empty()) continue;
    v_inputFileList.push_back(line);
  }
  inputFile.close();
  if (v_inputFileList.size() == 0) {
    std::cout << "No inputs provided" << std::endl;
    exit(1);
  }

  return;
}

void RecoTrainingNtupleMaker::GetBranches() {

  delphes_chain->SetBranchAddress("EventNumber", &evt);
  delphes_chain->SetBranchAddress("PassSelection", &passSelection);
  delphes_chain->SetBranchAddress("NumberOfJets", &numJets);
  delphes_chain->SetBranchAddress("NumberBTaggedJets", &numBTags);
  delphes_chain->SetBranchAddress("NumberMatchedJets", &numMatched);

  delphes_chain->SetBranchAddress("HasMatchSIG_bFromTopAndHiggs_bFromTop_bFromHiggs", &hasRequiredSIGMatches);
  delphes_chain->SetBranchAddress("HasMatchBKG_bFromTop1_bFromTop2", &hasRequiredBKGMatches);

  //-------------
  // Reco information
  //-------------
  v_AllJets_pt   = new std::vector<float>();
  v_AllJets_eta  = new std::vector<float>();
  v_AllJets_phi  = new std::vector<float>();
  v_AllJets_E    = new std::vector<float>();
  v_AllJets_M    = new std::vector<float>();
  v_AllJets_bTag = new std::vector<int>();
  delphes_chain->SetBranchAddress("Delphes_AllJets_pt", &v_AllJets_pt);
  delphes_chain->SetBranchAddress("Delphes_AllJets_eta", &v_AllJets_eta);
  delphes_chain->SetBranchAddress("Delphes_AllJets_phi", &v_AllJets_phi);
  delphes_chain->SetBranchAddress("Delphes_AllJets_E", &v_AllJets_E);
  delphes_chain->SetBranchAddress("Delphes_AllJets_M", &v_AllJets_M);
  delphes_chain->SetBranchAddress("Delphes_AllJets_BTag", &v_AllJets_bTag);

  for (int i = 0; i < 4; i++) v_SelectedJets_Lorentz.push_back(new TLorentzVector());
  delphes_chain->SetBranchAddress("Delphes_SelectedJet0_Lorentz", &v_SelectedJets_Lorentz.at(0));
  delphes_chain->SetBranchAddress("Delphes_SelectedJet1_Lorentz", &v_SelectedJets_Lorentz.at(1));
  delphes_chain->SetBranchAddress("Delphes_SelectedJet2_Lorentz", &v_SelectedJets_Lorentz.at(2));
  delphes_chain->SetBranchAddress("Delphes_SelectedJet3_Lorentz", &v_SelectedJets_Lorentz.at(3));

  v_SelectedJets_pt   = new std::vector<float>();
  v_SelectedJets_eta  = new std::vector<float>();
  v_SelectedJets_phi  = new std::vector<float>();
  v_SelectedJets_E    = new std::vector<float>();
  v_SelectedJets_M    = new std::vector<float>();
  v_SelectedJets_bTag = new std::vector<int>();
  delphes_chain->SetBranchAddress("Delphes_SelectedJets_pt", &v_SelectedJets_pt);
  delphes_chain->SetBranchAddress("Delphes_SelectedJets_eta", &v_SelectedJets_eta);
  delphes_chain->SetBranchAddress("Delphes_SelectedJets_phi", &v_SelectedJets_phi);
  delphes_chain->SetBranchAddress("Delphes_SelectedJets_E", &v_SelectedJets_E);
  delphes_chain->SetBranchAddress("Delphes_SelectedJets_M", &v_SelectedJets_M);
  delphes_chain->SetBranchAddress("Delphes_SelectedJets_BTag", &v_SelectedJets_bTag);

  v_bIndexMatchedToJet = new std::vector<int>();
  v_JetMatchedOrigin = new std::vector<std::string>();
  delphes_chain->SetBranchAddress("MatchingInfo_IndexOfQuarkMatchedToJet", &v_bIndexMatchedToJet);
  delphes_chain->SetBranchAddress("MatchingInfo_OriginOfQuarkMatchedToJet", &v_JetMatchedOrigin);

  delphes_poslep = new TLorentzVector();
  delphes_chain->SetBranchAddress("Delphes_PosLep_Lorentz", &delphes_poslep);
  delphes_chain->SetBranchAddress("Delphes_PosLep_pt", &delphes_poslep_pt);
  delphes_chain->SetBranchAddress("Delphes_PosLep_eta", &delphes_poslep_eta);
  delphes_chain->SetBranchAddress("Delphes_PosLep_phi", &delphes_poslep_phi);
  delphes_chain->SetBranchAddress("Delphes_PosLep_E", &delphes_poslep_E);
  delphes_chain->SetBranchAddress("Delphes_PosLep_M", &delphes_poslep_M);

  delphes_neglep = new TLorentzVector();
  delphes_chain->SetBranchAddress("Delphes_NegLep_Lorentz", &delphes_neglep);
  delphes_chain->SetBranchAddress("Delphes_NegLep_pt", &delphes_neglep_pt);
  delphes_chain->SetBranchAddress("Delphes_NegLep_eta", &delphes_neglep_eta);
  delphes_chain->SetBranchAddress("Delphes_NegLep_phi", &delphes_neglep_phi);
  delphes_chain->SetBranchAddress("Delphes_NegLep_E", &delphes_neglep_E);
  delphes_chain->SetBranchAddress("Delphes_NegLep_M", &delphes_neglep_M);

  delphes_etmiss = new TLorentzVector();
  delphes_chain->SetBranchAddress("Delphes_ETmiss_Lorentz", &delphes_etmiss);
  delphes_chain->SetBranchAddress("Delphes_ETmiss_MET", &delphes_etmiss_met);
  delphes_chain->SetBranchAddress("Delphes_ETmiss_eta", &delphes_etmiss_eta);
  delphes_chain->SetBranchAddress("Delphes_ETmiss_phi", &delphes_etmiss_phi);

}

void RecoTrainingNtupleMaker::Run() {

  // Skimmed delphes ntuple
  delphes_chain = new TChain("RecoResults");
  for (auto file : v_inputFileList) {
    delphes_chain->AddFile(file.c_str());
  }
  // get objects from skimmed ntuple
  GetBranches();

  // training ntuple
  TFile* outfile = TFile::Open((m_outputFileName).c_str(), "RECREATE");
  outfile->cd();
  TTree* training_tree = new TTree("RecoResults", "RecoResults");
  // identifying event
  training_tree->Branch("EventNumber", &m_evtNum, "EventNumber/I");
  training_tree->Branch("NumberOfJets", &m_nJets, "NumberOfJets/I");
  training_tree->Branch("NumberBTaggedJets", &m_nBTags, "NumberBTaggedJets/I");
  // identifying permutation
  training_tree->Branch("Index_HcJet", &m_HcJet_index, "Index_HcJet/I");
  training_tree->Branch("Index_tHcJet", &m_tHcJet_index, "Index_tHcJet/I");
  training_tree->Branch("Index_tJet", &m_tJet_index, "Index_tJet/I");
  training_tree->Branch("Index_gJet", &m_gJet_index, "Index_gJet/I");
  // training variables
  training_tree->Branch("CorrectMatch", &m_correctmatch, "CorrectMatch/I");
  training_tree->Branch("Pt_HcJet", &m_HcJet_pt, "Pt_HcJet/F");
  training_tree->Branch("Pt_tHcJet", &m_tHcJet_pt, "Pt_tHcJet/F");
  training_tree->Branch("Pt_tJet", &m_tJet_pt, "Pt_tJet/F");
  training_tree->Branch("Pt_gJet", &m_gJet_pt, "Pt_gJet/F");
  training_tree->Branch("Eta_HcJet", &m_HcJet_eta, "Eta_HcJet/F");
  training_tree->Branch("Eta_tHcJet", &m_tHcJet_eta, "Eta_tHcJet/F");
  training_tree->Branch("Eta_tJet", &m_tJet_eta, "Eta_tJet/F");
  training_tree->Branch("Eta_gJet", &m_gJet_eta, "Eta_gJet/F");
  training_tree->Branch("Phi_HcJet", &m_HcJet_phi, "Phi_HcJet/F");
  training_tree->Branch("Phi_tHcJet", &m_tHcJet_phi, "Phi_tHcJet/F");
  training_tree->Branch("Phi_tJet", &m_tJet_phi, "Phi_tJet/F");
  training_tree->Branch("Phi_gJet", &m_gJet_phi, "Phi_gJet/F");
  training_tree->Branch("E_HcJet", &m_HcJet_e, "E_HcJet/F");
  training_tree->Branch("E_tHcJet", &m_tHcJet_e, "E_tHcJet/F");
  training_tree->Branch("E_tJet", &m_tJet_e, "E_tJet/F");
  training_tree->Branch("E_gJet", &m_gJet_e, "E_gJet/F");
  training_tree->Branch("dPtPair_HcJet_tHcJet", &m_HcJet_tHcJet_dpt, "dPtPair_HcJet_tHcJet/F");
  training_tree->Branch("dPtPair_HcJet_tJet", &m_HcJet_tJet_dpt, "dPtPair_HcJet_tJet/F");
  training_tree->Branch("dPtPair_HcJet_gJet", &m_HcJet_gJet_dpt, "dPtPair_HcJet_gJet/F");
  training_tree->Branch("dPtPair_tHcJet_tJet", &m_tHcJet_tJet_dpt, "dPtPair_tHcJet_tJet/F");
  training_tree->Branch("dPtPair_tHcJet_gJet", &m_tHcJet_gJet_dpt, "dPtPair_tHcJet_gJet/F");
  training_tree->Branch("dPtPair_tJet_gJet", &m_tJet_gJet_dpt, "dPtPair_tJet_gJet/F");
  training_tree->Branch("dEtaPair_HcJet_tHcJet", &m_HcJet_tHcJet_deta, "dEtaPair_HcJet_tHcJet/F");
  training_tree->Branch("dEtaPair_HcJet_tJet", &m_HcJet_tJet_deta, "dEtaPair_HcJet_tJet/F");
  training_tree->Branch("dEtaPair_HcJet_gJet", &m_HcJet_gJet_deta, "dEtaPair_HcJet_gJet/F");
  training_tree->Branch("dEtaPair_tHcJet_tJet", &m_tHcJet_tJet_deta, "dEtaPair_tHcJet_tJet/F");
  training_tree->Branch("dEtaPair_tHcJet_gJet", &m_tHcJet_gJet_deta, "dEtaPair_tHcJet_gJet/F");
  training_tree->Branch("dEtaPair_tJet_gJet", &m_tJet_gJet_deta, "dEtaPair_tJet_gJet/F");
  training_tree->Branch("dPhiPair_HcJet_tHcJet", &m_HcJet_tHcJet_dphi, "dPhiPair_HcJet_tHcJet/F");
  training_tree->Branch("dPhiPair_HcJet_tJet", &m_HcJet_tJet_dphi, "dPhiPair_HcJet_tJet/F");
  training_tree->Branch("dPhiPair_HcJet_gJet", &m_HcJet_gJet_dphi, "dPhiPair_HcJet_gJet/F");
  training_tree->Branch("dPhiPair_tHcJet_tJet", &m_tHcJet_tJet_dphi, "dPhiPair_tHcJet_tJet/F");
  training_tree->Branch("dPhiPair_tHcJet_gJet", &m_tHcJet_gJet_dphi, "dPhiPair_tHcJet_gJet/F");
  training_tree->Branch("dPhiPair_tJet_gJet", &m_tJet_gJet_dphi, "dPhiPair_tJet_gJet/F");
  training_tree->Branch("dMPair_HcJet_tHcJet", &m_HcJet_tHcJet_dm, "dMPair_HcJet_tHcJet/F");
  training_tree->Branch("dMPair_HcJet_tJet", &m_HcJet_tJet_dm, "dMPair_HcJet_tJet/F");
  training_tree->Branch("dMPair_HcJet_gJet", &m_HcJet_gJet_dm, "dMPair_HcJet_gJet/F");
  training_tree->Branch("dMPair_tHcJet_tJet", &m_tHcJet_tJet_dm, "dMPair_tHcJet_tJet/F");
  training_tree->Branch("dMPair_tHcJet_gJet", &m_tHcJet_gJet_dm, "dMPair_tHcJet_gJet/F");
  training_tree->Branch("dMPair_tJet_gJet", &m_tJet_gJet_dm, "dMPair_tJet_gJet/F");
  training_tree->Branch("dRPair_HcJet_tHcJet", &m_HcJet_tHcJet_dR, "dRPair_HcJet_tHcJet/F");
  training_tree->Branch("dRPair_HcJet_tJet", &m_HcJet_tJet_dR, "dRPair_HcJet_tJet/F");
  training_tree->Branch("dRPair_HcJet_gJet", &m_HcJet_gJet_dR, "dRPair_HcJet_gJet/F");
  training_tree->Branch("dRPair_tHcJet_tJet", &m_tHcJet_tJet_dR, "dRPair_tHcJet_tJet/F");
  training_tree->Branch("dRPair_tHcJet_gJet", &m_tHcJet_gJet_dR, "dRPair_tHcJet_gJet/F");
  training_tree->Branch("dRPair_tJet_gJet", &m_tJet_gJet_dR, "dRPair_tJet_gJet/F");
  training_tree->Branch("PtPair_HcJet_tHcJet", &m_HcJet_tHcJet_ptpair, "PtPair_HcJet_tHcJet/F");
  training_tree->Branch("PtPair_HcJet_tJet", &m_HcJet_tJet_ptpair, "PtPair_HcJet_tJet/F");
  training_tree->Branch("PtPair_HcJet_gJet", &m_HcJet_gJet_ptpair, "PtPair_HcJet_gJet/F");
  training_tree->Branch("PtPair_tHcJet_tJet", &m_tHcJet_tJet_ptpair, "PtPair_tHcJet_tJet/F");
  training_tree->Branch("PtPair_tHcJet_gJet", &m_tHcJet_gJet_ptpair, "PtPair_tHcJet_gJet/F");
  training_tree->Branch("PtPair_tJet_gJet", &m_tJet_gJet_ptpair, "PtPair_tJet_gJet/F");
  training_tree->Branch("EtaPair_HcJet_tHcJet", &m_HcJet_tHcJet_etapair, "EtaPair_HcJet_tHcJet/F");
  training_tree->Branch("EtaPair_HcJet_tJet", &m_HcJet_tJet_etapair, "EtaPair_HcJet_tJet/F");
  training_tree->Branch("EtaPair_HcJet_gJet", &m_HcJet_gJet_etapair, "EtaPair_HcJet_gJet/F");
  training_tree->Branch("EtaPair_tHcJet_tJet", &m_tHcJet_tJet_etapair, "EtaPair_tHcJet_tJet/F");
  training_tree->Branch("EtaPair_tHcJet_gJet", &m_tHcJet_gJet_etapair, "EtaPair_tHcJet_gJet/F");
  training_tree->Branch("EtaPair_tJet_gJet", &m_tJet_gJet_etapair, "EtaPair_tJet_gJet/F");
  training_tree->Branch("PhiPair_HcJet_tHcJet", &m_HcJet_tHcJet_phipair, "PhiPair_HcJet_tHcJet/F");
  training_tree->Branch("PhiPair_HcJet_tJet", &m_HcJet_tJet_phipair, "PhiPair_HcJet_tJet/F");
  training_tree->Branch("PhiPair_HcJet_gJet", &m_HcJet_gJet_phipair, "PhiPair_HcJet_gJet/F");
  training_tree->Branch("PhiPair_tHcJet_tJet", &m_tHcJet_tJet_phipair, "PhiPair_tHcJet_tJet/F");
  training_tree->Branch("PhiPair_tHcJet_gJet", &m_tHcJet_gJet_phipair, "PhiPair_tHcJet_gJet/F");
  training_tree->Branch("PhiPair_tJet_gJet", &m_tJet_gJet_phipair, "PhiPair_tJet_gJet/F");
  training_tree->Branch("EPair_HcJet_tHcJet", &m_HcJet_tHcJet_epair, "EPair_HcJet_tHcJet/F");
  training_tree->Branch("EPair_HcJet_tJet", &m_HcJet_tJet_epair, "EPair_HcJet_tJet/F");
  training_tree->Branch("EPair_HcJet_gJet", &m_HcJet_gJet_epair, "EPair_HcJet_gJet/F");
  training_tree->Branch("EPair_tHcJet_tJet", &m_tHcJet_tJet_epair, "EPair_tHcJet_tJet/F");
  training_tree->Branch("EPair_tHcJet_gJet", &m_tHcJet_gJet_epair, "EPair_tHcJet_gJet/F");
  training_tree->Branch("EPair_tJet_gJet", &m_tJet_gJet_epair, "EPair_tJet_gJet/F");
  training_tree->Branch("Minv_HcJet_PosLep", &m_HcJet_PosLep_minv, "Minv_HcJet_PosLep/F");
  training_tree->Branch("Minv_tHcJet_PosLep", &m_tHcJet_PosLep_minv, "Minv_tHcJet_PosLep/F");
  training_tree->Branch("Minv_tJet_PosLep", &m_tJet_PosLep_minv, "Minv_tJet_PosLep/F");
  training_tree->Branch("Minv_gJet_PosLep", &m_gJet_PosLep_minv, "Minv_gJet_PosLep/F");
  training_tree->Branch("Minv_HcJet_NegLep", &m_HcJet_NegLep_minv, "Minv_HcJet_NegLep/F");
  training_tree->Branch("Minv_tHcJet_NegLep", &m_tHcJet_NegLep_minv, "Minv_tHcJet_NegLep/F");
  training_tree->Branch("Minv_tJet_NegLep", &m_tJet_NegLep_minv, "Minv_tJet_NegLep/F");
  training_tree->Branch("Minv_gJet_NegLep", &m_gJet_NegLep_minv, "Minv_gJet_NegLep/F");
  training_tree->Branch("PtPair_HcJet_PosLep", &m_HcJet_PosLep_ptpair, "PtPair_HcJet_PosLep/F");
  training_tree->Branch("PtPair_tHcJet_PosLep", &m_tHcJet_PosLep_ptpair, "PtPair_tHcJet_PosLep/F");
  training_tree->Branch("PtPair_tJet_PosLep", &m_tJet_PosLep_ptpair, "PtPair_tJet_PosLep/F");
  training_tree->Branch("PtPair_gJet_PosLep", &m_gJet_PosLep_ptpair, "PtPair_gJet_PosLep/F");
  training_tree->Branch("PtPair_HcJet_NegLep", &m_HcJet_NegLep_ptpair, "PtPair_HcJet_NegLep/F");
  training_tree->Branch("PtPair_tHcJet_NegLep", &m_tHcJet_NegLep_ptpair, "PtPair_tHcJet_NegLep/F");
  training_tree->Branch("PtPair_tJet_NegLep", &m_tJet_NegLep_ptpair, "PtPair_tJet_NegLep/F");
  training_tree->Branch("PtPair_gJet_NegLep", &m_gJet_NegLep_ptpair, "PtPair_gJet_NegLep/F");
  training_tree->Branch("EtaPair_HcJet_PosLep", &m_HcJet_PosLep_etapair, "EtaPair_HcJet_PosLep/F");
  training_tree->Branch("EtaPair_tHcJet_PosLep", &m_tHcJet_PosLep_etapair, "EtaPair_tHcJet_PosLep/F");
  training_tree->Branch("EtaPair_tJet_PosLep", &m_tJet_PosLep_etapair, "EtaPair_tJet_PosLep/F");
  training_tree->Branch("EtaPair_gJet_PosLep", &m_gJet_PosLep_etapair, "EtaPair_gJet_PosLep/F");
  training_tree->Branch("EtaPair_HcJet_NegLep", &m_HcJet_NegLep_etapair, "EtaPair_HcJet_NegLep/F");
  training_tree->Branch("EtaPair_tHcJet_NegLep", &m_tHcJet_NegLep_etapair, "EtaPair_tHcJet_NegLep/F");
  training_tree->Branch("EtaPair_tJet_NegLep", &m_tJet_NegLep_etapair, "EtaPair_tJet_NegLep/F");
  training_tree->Branch("EtaPair_gJet_NegLep", &m_gJet_NegLep_etapair, "EtaPair_gJet_NegLep/F");
  training_tree->Branch("PhiPair_HcJet_PosLep", &m_HcJet_PosLep_phipair, "PhiPair_HcJet_PosLep/F");
  training_tree->Branch("PhiPair_tHcJet_PosLep", &m_tHcJet_PosLep_phipair, "PhiPair_tHcJet_PosLep/F");
  training_tree->Branch("PhiPair_tJet_PosLep", &m_tJet_PosLep_phipair, "PhiPair_tJet_PosLep/F");
  training_tree->Branch("PhiPair_gJet_PosLep", &m_gJet_PosLep_phipair, "PhiPair_gJet_PosLep/F");
  training_tree->Branch("PhiPair_HcJet_NegLep", &m_HcJet_NegLep_phipair, "PhiPair_HcJet_NegLep/F");
  training_tree->Branch("PhiPair_tHcJet_NegLep", &m_tHcJet_NegLep_phipair, "PhiPair_tHcJet_NegLep/F");
  training_tree->Branch("PhiPair_tJet_NegLep", &m_tJet_NegLep_phipair, "PhiPair_tJet_NegLep/F");
  training_tree->Branch("PhiPair_gJet_NegLep", &m_gJet_NegLep_phipair, "PhiPair_gJet_NegLep/F");
  training_tree->Branch("EPair_HcJet_PosLep", &m_HcJet_PosLep_epair, "EPair_HcJet_PosLep/F");
  training_tree->Branch("EPair_tHcJet_PosLep", &m_tHcJet_PosLep_epair, "EPair_tHcJet_PosLep/F");
  training_tree->Branch("EPair_tJet_PosLep", &m_tJet_PosLep_epair, "EPair_tJet_PosLep/F");
  training_tree->Branch("EPair_gJet_PosLep", &m_gJet_PosLep_epair, "EPair_gJet_PosLep/F");
  training_tree->Branch("EPair_HcJet_NegLep", &m_HcJet_NegLep_epair, "EPair_HcJet_NegLep/F");
  training_tree->Branch("EPair_tHcJet_NegLep", &m_tHcJet_NegLep_epair, "EPair_tHcJet_NegLep/F");
  training_tree->Branch("EPair_tJet_NegLep", &m_tJet_NegLep_epair, "EPair_tJet_NegLep/F");
  training_tree->Branch("EPair_gJet_NegLep", &m_gJet_NegLep_epair, "EPair_gJet_NegLep/F");

  // run over events
  int nEntries = delphes_chain->GetEntries();
  for (int entry = 0; entry < nEntries; entry++) {
    delphes_chain->GetEntry(entry);
    m_evtNum = entry;
    m_nJets  = numJets;
    m_nBTags = numBTags;

    if (entry%10000 == 0) std::cout << "Event " << entry << std::endl;

    if (passSelection == 1) {
      // order is jet taken as: Hc, topHc, top, gluon
      for (int Hc = 0; Hc < 4; ++Hc) {
        for (int tHc = 0; tHc < 4; ++tHc) {
          if (tHc == Hc) continue;
          for (int t = 0; t < 4; ++t) {
            if (t == Hc || t == tHc) continue;
            for (int g = 0; g < 4; ++g) {
              if (g == Hc || g == tHc || g == t) continue;

              m_correctmatch = 0;

              // check if Hc is bFromHiggs AND tHc is bFromTopAndHiggs AND t is bFromTop
              if (v_JetMatchedOrigin->at(Hc).compare("bFromHiggs") == 0 && v_JetMatchedOrigin->at(tHc).compare("bFromTopAndHiggs") == 0 && v_JetMatchedOrigin->at(t).compare("bFromTop") == 0) {
                m_correctmatch = 1;
              }

              m_HcJet_index  = Hc;
              m_tHcJet_index = tHc;
              m_tJet_index   = t;
              m_gJet_index   = g;

              TLorentzVector* l_Hc  = v_SelectedJets_Lorentz.at(Hc);
              TLorentzVector* l_tHc = v_SelectedJets_Lorentz.at(tHc);
              TLorentzVector* l_top = v_SelectedJets_Lorentz.at(t);
              TLorentzVector* l_g   = v_SelectedJets_Lorentz.at(g);

              m_HcJet_pt  = l_Hc->Pt();
              m_tHcJet_pt = l_tHc->Pt();
              m_tJet_pt   = l_top->Pt();
              m_gJet_pt   = l_g->Pt();

              m_HcJet_eta  = l_Hc->Eta();
              m_tHcJet_eta = l_tHc->Eta();
              m_tJet_eta   = l_top->Eta();
              m_gJet_eta   = l_g->Eta();

              m_HcJet_phi  = l_Hc->Phi();
              m_tHcJet_phi = l_tHc->Phi();
              m_tJet_phi   = l_top->Phi();
              m_gJet_phi   = l_g->Phi();

              m_HcJet_e  = l_Hc->E();
              m_tHcJet_e = l_tHc->E();
              m_tJet_e   = l_top->E();
              m_gJet_e   = l_g->E();

              m_HcJet_tHcJet_dpt = DPtPair(l_Hc, l_tHc);
              m_HcJet_tJet_dpt   = DPtPair(l_Hc, l_top);
              m_HcJet_gJet_dpt   = DPtPair(l_Hc, l_g);
              m_tHcJet_tJet_dpt  = DPtPair(l_tHc, l_top);
              m_tHcJet_gJet_dpt  = DPtPair(l_tHc, l_g);
              m_tJet_gJet_dpt    = DPtPair(l_top, l_g);

              m_HcJet_tHcJet_deta = DEtaPair(l_Hc, l_tHc);
              m_HcJet_tJet_deta   = DEtaPair(l_Hc, l_top);
              m_HcJet_gJet_deta   = DEtaPair(l_Hc, l_g);
              m_tHcJet_tJet_deta  = DEtaPair(l_tHc, l_top);
              m_tHcJet_gJet_deta  = DEtaPair(l_tHc, l_g);
              m_tJet_gJet_deta    = DEtaPair(l_top, l_g);

              m_HcJet_tHcJet_dphi = DPhiPair(l_Hc, l_tHc);
              m_HcJet_tJet_dphi   = DPhiPair(l_Hc, l_top);
              m_HcJet_gJet_dphi   = DPhiPair(l_Hc, l_g);
              m_tHcJet_tJet_dphi  = DPhiPair(l_tHc, l_top);
              m_tHcJet_gJet_dphi  = DPhiPair(l_tHc, l_g);
              m_tJet_gJet_dphi    = DPhiPair(l_top, l_g);

              m_HcJet_tHcJet_dm = DMPair(l_Hc, l_tHc);
              m_HcJet_tJet_dm   = DMPair(l_Hc, l_top);
              m_HcJet_gJet_dm   = DMPair(l_Hc, l_g);
              m_tHcJet_tJet_dm  = DMPair(l_tHc, l_top);
              m_tHcJet_gJet_dm  = DMPair(l_tHc, l_g);
              m_tJet_gJet_dm    = DMPair(l_top, l_g);

              m_HcJet_tHcJet_dR = DRPair(l_Hc, l_tHc);
              m_HcJet_tJet_dR   = DRPair(l_Hc, l_top);
              m_HcJet_gJet_dR   = DRPair(l_Hc, l_g);
              m_tHcJet_tJet_dR  = DRPair(l_tHc, l_top);
              m_tHcJet_gJet_dR  = DRPair(l_tHc, l_g);
              m_tJet_gJet_dR    = DRPair(l_top, l_g);

              m_HcJet_tHcJet_ptpair = PtPair(l_Hc, l_tHc);
              m_HcJet_tJet_ptpair   = PtPair(l_Hc, l_top);
              m_HcJet_gJet_ptpair   = PtPair(l_Hc, l_g);
              m_tHcJet_tJet_ptpair  = PtPair(l_tHc, l_top);
              m_tHcJet_gJet_ptpair  = PtPair(l_tHc, l_g);
              m_tJet_gJet_ptpair    = PtPair(l_top, l_g);

              m_HcJet_tHcJet_etapair = EtaPair(l_Hc, l_tHc);
              m_HcJet_tJet_etapair   = EtaPair(l_Hc, l_top);
              m_HcJet_gJet_etapair   = EtaPair(l_Hc, l_g);
              m_tHcJet_tJet_etapair  = EtaPair(l_tHc, l_top);
              m_tHcJet_gJet_etapair  = EtaPair(l_tHc, l_g);
              m_tJet_gJet_etapair    = EtaPair(l_top, l_g);

              m_HcJet_tHcJet_phipair = PhiPair(l_Hc, l_tHc);
              m_HcJet_tJet_phipair   = PhiPair(l_Hc, l_top);
              m_HcJet_gJet_phipair   = PhiPair(l_Hc, l_g);
              m_tHcJet_tJet_phipair  = PhiPair(l_tHc, l_top);
              m_tHcJet_gJet_phipair  = PhiPair(l_tHc, l_g);
              m_tJet_gJet_phipair    = PhiPair(l_top, l_g);

              m_HcJet_tHcJet_epair = EPair(l_Hc, l_tHc);
              m_HcJet_tJet_epair   = EPair(l_Hc, l_top);
              m_HcJet_gJet_epair   = EPair(l_Hc, l_g);
              m_tHcJet_tJet_epair  = EPair(l_tHc, l_top);
              m_tHcJet_gJet_epair  = EPair(l_tHc, l_g);
              m_tJet_gJet_epair    = EPair(l_top, l_g);

              m_HcJet_PosLep_minv  = MinvPair(l_Hc, delphes_poslep);
              m_tHcJet_PosLep_minv = MinvPair(l_tHc, delphes_poslep);
              m_tJet_PosLep_minv   = MinvPair(l_top, delphes_poslep);
              m_gJet_PosLep_minv   = MinvPair(l_g, delphes_poslep);
              m_HcJet_NegLep_minv  = MinvPair(l_Hc, delphes_neglep);
              m_tHcJet_NegLep_minv = MinvPair(l_tHc, delphes_neglep);
              m_tJet_NegLep_minv   = MinvPair(l_top, delphes_neglep);
              m_gJet_NegLep_minv   = MinvPair(l_g, delphes_neglep);

              m_HcJet_PosLep_ptpair  = PtPair(l_Hc, delphes_poslep);
              m_tHcJet_PosLep_ptpair = PtPair(l_tHc, delphes_poslep);
              m_tJet_PosLep_ptpair   = PtPair(l_top, delphes_poslep);
              m_gJet_PosLep_ptpair   = PtPair(l_g, delphes_poslep);
              m_HcJet_NegLep_ptpair  = PtPair(l_Hc, delphes_neglep);
              m_tHcJet_NegLep_ptpair = PtPair(l_tHc, delphes_neglep);
              m_tJet_NegLep_ptpair   = PtPair(l_top, delphes_neglep);
              m_gJet_NegLep_ptpair   = PtPair(l_g, delphes_neglep);

              m_HcJet_PosLep_etapair  = EtaPair(l_Hc, delphes_poslep);
              m_tHcJet_PosLep_etapair = EtaPair(l_tHc, delphes_poslep);
              m_tJet_PosLep_etapair   = EtaPair(l_top, delphes_poslep);
              m_gJet_PosLep_etapair   = EtaPair(l_g, delphes_poslep);
              m_HcJet_NegLep_etapair  = EtaPair(l_Hc, delphes_neglep);
              m_tHcJet_NegLep_etapair = EtaPair(l_tHc, delphes_neglep);
              m_tJet_NegLep_etapair   = EtaPair(l_top, delphes_neglep);
              m_gJet_NegLep_etapair   = EtaPair(l_g, delphes_neglep);

              m_HcJet_PosLep_phipair  = PhiPair(l_Hc, delphes_poslep);
              m_tHcJet_PosLep_phipair = PhiPair(l_tHc, delphes_poslep);
              m_tJet_PosLep_phipair   = PhiPair(l_top, delphes_poslep);
              m_gJet_PosLep_phipair   = PhiPair(l_g, delphes_poslep);
              m_HcJet_NegLep_phipair  = PhiPair(l_Hc, delphes_neglep);
              m_tHcJet_NegLep_phipair = PhiPair(l_tHc, delphes_neglep);
              m_tJet_NegLep_phipair   = PhiPair(l_top, delphes_neglep);
              m_gJet_NegLep_phipair   = PhiPair(l_g, delphes_neglep);

              m_HcJet_PosLep_epair  = EPair(l_Hc, delphes_poslep);
              m_tHcJet_PosLep_epair = EPair(l_tHc, delphes_poslep);
              m_tJet_PosLep_epair   = EPair(l_top, delphes_poslep);
              m_gJet_PosLep_epair   = EPair(l_g, delphes_poslep);
              m_HcJet_NegLep_epair  = EPair(l_Hc, delphes_neglep);
              m_tHcJet_NegLep_epair = EPair(l_tHc, delphes_neglep);
              m_tJet_NegLep_epair   = EPair(l_top, delphes_neglep);
              m_gJet_NegLep_epair   = EPair(l_g, delphes_neglep);

              training_tree->Fill();
            }
          }
        }
      }
    }
  }

  training_tree->Write("", TObject::kOverwrite);
  delete training_tree;
  outfile->Close();

  delete delphes_chain;
  DeleteChainVariables();

  return;
}

float RecoTrainingNtupleMaker::DPtPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->Pt() - b->Pt());
}

float RecoTrainingNtupleMaker::DEtaPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->Eta() - b->Eta());
}

float RecoTrainingNtupleMaker::DPhiPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->DeltaPhi(*b));
}

float RecoTrainingNtupleMaker::DMPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->M() - b->M());
}

float RecoTrainingNtupleMaker::DRPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->DeltaR(*b));
}

float RecoTrainingNtupleMaker::PtPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).Pt();
}

float RecoTrainingNtupleMaker::MinvPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).M();
}

float RecoTrainingNtupleMaker::EtaPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).Eta();
}

float RecoTrainingNtupleMaker::PhiPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).Phi();
}

float RecoTrainingNtupleMaker::EPair(TLorentzVector* a, TLorentzVector* b) {
   return ((*a)+(*b)).E();
}

void RecoTrainingNtupleMaker::DeleteChainVariables() {

  for (auto v : v_SelectedJets_Lorentz) delete v;

  v_SelectedJets_Lorentz.clear();
  v_SelectedJets_pt->clear();
  v_SelectedJets_eta->clear();
  v_SelectedJets_phi->clear();
  v_SelectedJets_E->clear();
  v_SelectedJets_M->clear();
  v_SelectedJets_bTag->clear();
  v_JetMatchedOrigin->clear();
  v_bIndexMatchedToJet->clear();

  delete v_SelectedJets_pt;
  delete v_SelectedJets_eta;
  delete v_SelectedJets_phi;
  delete v_SelectedJets_E;
  delete v_SelectedJets_M;
  delete v_SelectedJets_bTag;
  delete v_JetMatchedOrigin;
  delete v_bIndexMatchedToJet;

  delete delphes_poslep;
  delete delphes_neglep;
  delete delphes_etmiss;
}
