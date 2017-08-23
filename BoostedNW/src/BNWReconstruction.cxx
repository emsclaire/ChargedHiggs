#include "include/BNWReconstruction.h"
#include "include/NeutrinoWeighter.h"

BNWReconstruction::BNWReconstruction(std::string inputFileName, std::string outputFileName, std::string typeBDT) :
m_inputFileName(inputFileName),
m_outputFileName(outputFileName) {

 m_WeightsDir    = "weights/";
 m_WeigthsPrefix = "TMVAClassification";
 m_typeBDT = typeBDT;

 m_WeightsFile = m_WeightsDir + m_WeigthsPrefix + "_BDT_"+typeBDT+".weights.xml";
 m_inputFile = new TFile(inputFileName.c_str());
 m_inputTree = (TTree*) m_inputFile->Get("RecoResults");
 m_outputFile = new TFile(outputFileName.c_str(),"RECREATE");

 BookNewTree();

}

void BNWReconstruction::InitialiseEventVariables() {

  maxNeutrinoWeight  = 0; // if stays 0, there is no solution / did not pass cuts
  isReconstructed    = 0; // 1 if NW finds solution
  MaxBDTWeight       = -99;


  maxNW_top     = TLorentzVector();
  maxNW_top_pt  = -99.;
  maxNW_top_eta = -99.;
  maxNW_top_phi = -99.;
  maxNW_top_E   = -99.;
  maxNW_top_M   = -99.;

  maxNW_tbar     = TLorentzVector();
  maxNW_tbar_pt  = -99.;
  maxNW_tbar_eta = -99.;
  maxNW_tbar_phi = -99.;
  maxNW_tbar_E   = -99.;
  maxNW_tbar_M   = -99.;

  maxNW_b     = TLorentzVector();
  maxNW_b_pt  = -99.;
  maxNW_b_eta = -99.;
  maxNW_b_phi = -99.;
  maxNW_b_E   = -99.;
  maxNW_b_M   = -99.;

  maxNW_bbar     = TLorentzVector();
  maxNW_bbar_pt  = -99.;
  maxNW_bbar_eta = -99.;
  maxNW_bbar_phi = -99.;
  maxNW_bbar_E   = -99.;
  maxNW_bbar_M   = -99.;



  maxNW_nu     = TLorentzVector();
  maxNW_nu_pt  = -99.;
  maxNW_nu_eta = -99.;
  maxNW_nu_phi = -99.;
  maxNW_nu_E   = -99.;
  maxNW_nu_M   = -99.;

  maxNW_nubar     = TLorentzVector();
  maxNW_nubar_pt  = -99.;
  maxNW_nubar_eta = -99.;
  maxNW_nubar_phi = -99.;
  maxNW_nubar_E   = -99.;
  maxNW_nubar_M   = -99.;



  maxNW_b_index    = -99;
  maxNW_bbar_index = -99;
  maxNW_bH_index = -99;


}


void BNWReconstruction::GetBranches() {

  m_inputTree->SetBranchAddress("EventNumber", &evt);
  m_inputTree->SetBranchAddress("PassSelection", &passSelection);
  m_inputTree->SetBranchAddress("NumberOfJets", &numJets);
  m_inputTree->SetBranchAddress("NumberBTaggedJets", &numBTags);
  m_inputTree->SetBranchAddress("NumberMatchedJets", &numMatched);

  m_inputTree->SetBranchAddress("HasMatchSIG_bFromTopAndHiggs_bFromTop_bFromHiggs", &hasRequiredSIGMatches);
  m_inputTree->SetBranchAddress("HasMatchBKG_bFromTop1_bFromTop2", &hasRequiredBKGMatches);

  //-------------
  // Reco information
  //-------------
  v_AllJets_pt   = new std::vector<float>();
  v_AllJets_eta  = new std::vector<float>();
  v_AllJets_phi  = new std::vector<float>();
  v_AllJets_E    = new std::vector<float>();
  v_AllJets_M    = new std::vector<float>();
  v_AllJets_bTag = new std::vector<int>();

  m_inputTree->SetBranchAddress("Delphes_AllJets_pt", &v_AllJets_pt);
  m_inputTree->SetBranchAddress("Delphes_AllJets_eta", &v_AllJets_eta);
  m_inputTree->SetBranchAddress("Delphes_AllJets_phi", &v_AllJets_phi);
  m_inputTree->SetBranchAddress("Delphes_AllJets_E", &v_AllJets_E);
  m_inputTree->SetBranchAddress("Delphes_AllJets_M", &v_AllJets_M);
  m_inputTree->SetBranchAddress("Delphes_AllJets_BTag", &v_AllJets_bTag);

  for (int i = 0; i < 4; i++) v_SelectedJets_Lorentz.push_back(new TLorentzVector());

    m_inputTree->SetBranchAddress("Delphes_SelectedJet0_Lorentz", &v_SelectedJets_Lorentz.at(0));
  m_inputTree->SetBranchAddress("Delphes_SelectedJet1_Lorentz", &v_SelectedJets_Lorentz.at(1));
  m_inputTree->SetBranchAddress("Delphes_SelectedJet2_Lorentz", &v_SelectedJets_Lorentz.at(2));
  m_inputTree->SetBranchAddress("Delphes_SelectedJet3_Lorentz", &v_SelectedJets_Lorentz.at(3));

  v_SelectedJets_pt   = new std::vector<float>();
  v_SelectedJets_eta  = new std::vector<float>();
  v_SelectedJets_phi  = new std::vector<float>();
  v_SelectedJets_E    = new std::vector<float>();
  v_SelectedJets_M    = new std::vector<float>();
  v_SelectedJets_bTag = new std::vector<int>();

  m_inputTree->SetBranchAddress("Delphes_SelectedJets_pt", &v_SelectedJets_pt);
  m_inputTree->SetBranchAddress("Delphes_SelectedJets_eta", &v_SelectedJets_eta);
  m_inputTree->SetBranchAddress("Delphes_SelectedJets_phi", &v_SelectedJets_phi);
  m_inputTree->SetBranchAddress("Delphes_SelectedJets_E", &v_SelectedJets_E);
  m_inputTree->SetBranchAddress("Delphes_SelectedJets_M", &v_SelectedJets_M);
  m_inputTree->SetBranchAddress("Delphes_SelectedJets_BTag", &v_SelectedJets_bTag);

  v_bIndexMatchedToJet = new std::vector<int>();
  v_JetMatchedOrigin = new std::vector<std::string>();

  m_inputTree->SetBranchAddress("MatchingInfo_IndexOfQuarkMatchedToJet", &v_bIndexMatchedToJet);
  m_inputTree->SetBranchAddress("MatchingInfo_OriginOfQuarkMatchedToJet", &v_JetMatchedOrigin);

  delphes_poslep = new TLorentzVector();

  m_inputTree->SetBranchAddress("Delphes_PosLep_Lorentz", &delphes_poslep);
  m_inputTree->SetBranchAddress("Delphes_PosLep_pt", &delphes_poslep_pt);
  m_inputTree->SetBranchAddress("Delphes_PosLep_eta", &delphes_poslep_eta);
  m_inputTree->SetBranchAddress("Delphes_PosLep_phi", &delphes_poslep_phi);
  m_inputTree->SetBranchAddress("Delphes_PosLep_E", &delphes_poslep_E);
  m_inputTree->SetBranchAddress("Delphes_PosLep_M", &delphes_poslep_M);

  delphes_neglep = new TLorentzVector();

  m_inputTree->SetBranchAddress("Delphes_NegLep_Lorentz", &delphes_neglep);
  m_inputTree->SetBranchAddress("Delphes_NegLep_pt", &delphes_neglep_pt);
  m_inputTree->SetBranchAddress("Delphes_NegLep_eta", &delphes_neglep_eta);
  m_inputTree->SetBranchAddress("Delphes_NegLep_phi", &delphes_neglep_phi);
  m_inputTree->SetBranchAddress("Delphes_NegLep_E", &delphes_neglep_E);
  m_inputTree->SetBranchAddress("Delphes_NegLep_M", &delphes_neglep_M);

  delphes_etmiss = new TLorentzVector();

  m_inputTree->SetBranchAddress("Delphes_ETmiss_Lorentz", &delphes_etmiss);
  m_inputTree->SetBranchAddress("Delphes_ETmiss_MET", &delphes_etmiss_met);
  m_inputTree->SetBranchAddress("Delphes_ETmiss_eta", &delphes_etmiss_eta);
  m_inputTree->SetBranchAddress("Delphes_ETmiss_phi", &delphes_etmiss_phi);

}

void BNWReconstruction::SetTMVAReader() {

  m_reader = new TMVA::Reader( "!Color:!Silent" );

  // m_reader->AddSpectator("EventNumber", &m_evtNum);
  // m_reader->AddSpectator("NumberOfJets", &m_nJets);
  // m_reader->AddSpectator("NumberBTaggedJets", &m_nBTags);
  // m_reader->AddSpectator("Index_HcJet", &m_HcJet_index);
  // m_reader->AddSpectator("Index_tHcJet", &m_tHcJet_index);
  // m_reader->AddSpectator("Index_tJet", &m_tJet_index);
  // m_reader->AddSpectator("Index_gJet", &m_gJet_index);
  // m_reader->AddSpectator("CorrectMatch",&m_correctmatch);


  m_reader->AddVariable("Pt_HcJet", &m_HcJet_pt);
  m_reader->AddVariable("Pt_tHcJet", &m_tHcJet_pt);
  m_reader->AddVariable("Pt_tJet", &m_tJet_pt);
  m_reader->AddVariable("Pt_gJet", &m_gJet_pt);
  m_reader->AddVariable("Eta_HcJet", &m_HcJet_eta);
  m_reader->AddVariable("Eta_tHcJet", &m_tHcJet_eta);
  m_reader->AddVariable("Eta_tJet", &m_tJet_eta);
  m_reader->AddVariable("Eta_gJet", &m_gJet_eta);
  m_reader->AddVariable("Phi_HcJet", &m_HcJet_phi);
  m_reader->AddVariable("Phi_tHcJet", &m_tHcJet_phi);
  m_reader->AddVariable("Phi_tJet", &m_tJet_phi);
  m_reader->AddVariable("Phi_gJet", &m_gJet_phi);
  m_reader->AddVariable("E_HcJet", &m_HcJet_e);
  m_reader->AddVariable("E_tHcJet", &m_tHcJet_e);
  m_reader->AddVariable("E_tJet", &m_tJet_e);
  m_reader->AddVariable("E_gJet", &m_gJet_e);
  m_reader->AddVariable("dPtPair_HcJet_tHcJet", &m_HcJet_tHcJet_dpt);
  m_reader->AddVariable("dPtPair_HcJet_tJet", &m_HcJet_tJet_dpt);
  m_reader->AddVariable("dPtPair_HcJet_gJet", &m_HcJet_gJet_dpt);
  m_reader->AddVariable("dPtPair_tHcJet_tJet", &m_tHcJet_tJet_dpt);
  m_reader->AddVariable("dPtPair_tHcJet_gJet", &m_tHcJet_gJet_dpt);
  m_reader->AddVariable("dPtPair_tJet_gJet", &m_tJet_gJet_dpt);
  m_reader->AddVariable("dEtaPair_HcJet_tHcJet", &m_HcJet_tHcJet_deta);
  m_reader->AddVariable("dEtaPair_HcJet_tJet", &m_HcJet_tJet_deta);
  m_reader->AddVariable("dEtaPair_HcJet_gJet", &m_HcJet_gJet_deta);
  m_reader->AddVariable("dEtaPair_tHcJet_tJet", &m_tHcJet_tJet_deta);
  m_reader->AddVariable("dEtaPair_tHcJet_gJet", &m_tHcJet_gJet_deta);
  m_reader->AddVariable("dEtaPair_tJet_gJet", &m_tJet_gJet_deta);
  m_reader->AddVariable("dPhiPair_HcJet_tHcJet", &m_HcJet_tHcJet_dphi);
  m_reader->AddVariable("dPhiPair_HcJet_tJet", &m_HcJet_tJet_dphi);
  m_reader->AddVariable("dPhiPair_HcJet_gJet", &m_HcJet_gJet_dphi);
  m_reader->AddVariable("dPhiPair_tHcJet_tJet", &m_tHcJet_tJet_dphi);
  m_reader->AddVariable("dPhiPair_tHcJet_gJet", &m_tHcJet_gJet_dphi);
  m_reader->AddVariable("dPhiPair_tJet_gJet", &m_tJet_gJet_dphi);
  m_reader->AddVariable("dMPair_HcJet_tHcJet", &m_HcJet_tHcJet_dm);
  m_reader->AddVariable("dMPair_HcJet_tJet", &m_HcJet_tJet_dm);
  m_reader->AddVariable("dMPair_HcJet_gJet", &m_HcJet_gJet_dm);
  m_reader->AddVariable("dMPair_tHcJet_tJet", &m_tHcJet_tJet_dm);
  m_reader->AddVariable("dMPair_tHcJet_gJet", &m_tHcJet_gJet_dm);
  m_reader->AddVariable("dMPair_tJet_gJet", &m_tJet_gJet_dm);
  m_reader->AddVariable("dRPair_HcJet_tHcJet", &m_HcJet_tHcJet_dR);
  m_reader->AddVariable("dRPair_HcJet_tJet", &m_HcJet_tJet_dR);
  m_reader->AddVariable("dRPair_HcJet_gJet", &m_HcJet_gJet_dR);
  m_reader->AddVariable("dRPair_tHcJet_tJet", &m_tHcJet_tJet_dR);
  m_reader->AddVariable("dRPair_tHcJet_gJet", &m_tHcJet_gJet_dR);
  m_reader->AddVariable("dRPair_tJet_gJet", &m_tJet_gJet_dR);
  m_reader->AddVariable("PtPair_HcJet_tHcJet", &m_HcJet_tHcJet_ptpair);
  m_reader->AddVariable("PtPair_HcJet_tJet", &m_HcJet_tJet_ptpair);
  m_reader->AddVariable("PtPair_HcJet_gJet", &m_HcJet_gJet_ptpair);
  m_reader->AddVariable("PtPair_tHcJet_tJet", &m_tHcJet_tJet_ptpair);
  m_reader->AddVariable("PtPair_tHcJet_gJet", &m_tHcJet_gJet_ptpair);
  m_reader->AddVariable("PtPair_tJet_gJet", &m_tJet_gJet_ptpair);
  m_reader->AddVariable("EtaPair_HcJet_tHcJet", &m_HcJet_tHcJet_etapair);
  m_reader->AddVariable("EtaPair_HcJet_tJet", &m_HcJet_tJet_etapair);
  m_reader->AddVariable("EtaPair_HcJet_gJet", &m_HcJet_gJet_etapair);
  m_reader->AddVariable("EtaPair_tHcJet_tJet", &m_tHcJet_tJet_etapair);
  m_reader->AddVariable("EtaPair_tHcJet_gJet", &m_tHcJet_gJet_etapair);
  m_reader->AddVariable("EtaPair_tJet_gJet", &m_tJet_gJet_etapair);
  m_reader->AddVariable("PhiPair_HcJet_tHcJet", &m_HcJet_tHcJet_phipair);
  m_reader->AddVariable("PhiPair_HcJet_tJet", &m_HcJet_tJet_phipair);
  m_reader->AddVariable("PhiPair_HcJet_gJet", &m_HcJet_gJet_phipair);
  m_reader->AddVariable("PhiPair_tHcJet_tJet", &m_tHcJet_tJet_phipair);
  m_reader->AddVariable("PhiPair_tHcJet_gJet", &m_tHcJet_gJet_phipair);
  m_reader->AddVariable("PhiPair_tJet_gJet", &m_tJet_gJet_phipair);
  m_reader->AddVariable("EPair_HcJet_tHcJet", &m_HcJet_tHcJet_epair);
  m_reader->AddVariable("EPair_HcJet_tJet", &m_HcJet_tJet_epair);
  m_reader->AddVariable("EPair_HcJet_gJet", &m_HcJet_gJet_epair);
  m_reader->AddVariable("EPair_tHcJet_tJet", &m_tHcJet_tJet_epair);
  m_reader->AddVariable("EPair_tHcJet_gJet", &m_tHcJet_gJet_epair);
  m_reader->AddVariable("EPair_tJet_gJet", &m_tJet_gJet_epair);
  m_reader->AddVariable("Minv_HcJet_PosLep", &m_HcJet_PosLep_minv);
  m_reader->AddVariable("Minv_tHcJet_PosLep", &m_tHcJet_PosLep_minv);
  m_reader->AddVariable("Minv_tJet_PosLep", &m_tJet_PosLep_minv);
  m_reader->AddVariable("Minv_gJet_PosLep", &m_gJet_PosLep_minv);
  m_reader->AddVariable("Minv_HcJet_NegLep", &m_HcJet_NegLep_minv);
  m_reader->AddVariable("Minv_tHcJet_NegLep", &m_tHcJet_NegLep_minv);
  m_reader->AddVariable("Minv_tJet_NegLep", &m_tJet_NegLep_minv);
  m_reader->AddVariable("Minv_gJet_NegLep", &m_gJet_NegLep_minv);
  m_reader->AddVariable("PtPair_HcJet_PosLep", &m_HcJet_PosLep_ptpair);
  m_reader->AddVariable("PtPair_tHcJet_PosLep", &m_tHcJet_PosLep_ptpair);
  m_reader->AddVariable("PtPair_tJet_PosLep", &m_tJet_PosLep_ptpair);
  m_reader->AddVariable("PtPair_gJet_PosLep", &m_gJet_PosLep_ptpair);
  m_reader->AddVariable("PtPair_HcJet_NegLep", &m_HcJet_NegLep_ptpair);
  m_reader->AddVariable("PtPair_tHcJet_NegLep", &m_tHcJet_NegLep_ptpair);
  m_reader->AddVariable("PtPair_tJet_NegLep", &m_tJet_NegLep_ptpair);
  m_reader->AddVariable("PtPair_gJet_NegLep", &m_gJet_NegLep_ptpair);
  m_reader->AddVariable("EtaPair_HcJet_PosLep", &m_HcJet_PosLep_etapair);
  m_reader->AddVariable("EtaPair_tHcJet_PosLep", &m_tHcJet_PosLep_etapair);
  m_reader->AddVariable("EtaPair_tJet_PosLep", &m_tJet_PosLep_etapair);
  m_reader->AddVariable("EtaPair_gJet_PosLep", &m_gJet_PosLep_etapair);
  m_reader->AddVariable("EtaPair_HcJet_NegLep", &m_HcJet_NegLep_etapair);
  m_reader->AddVariable("EtaPair_tHcJet_NegLep", &m_tHcJet_NegLep_etapair);
  m_reader->AddVariable("EtaPair_tJet_NegLep", &m_tJet_NegLep_etapair);
  m_reader->AddVariable("EtaPair_gJet_NegLep", &m_gJet_NegLep_etapair);
  m_reader->AddVariable("PhiPair_HcJet_PosLep", &m_HcJet_PosLep_phipair);
  m_reader->AddVariable("PhiPair_tHcJet_PosLep", &m_tHcJet_PosLep_phipair);
  m_reader->AddVariable("PhiPair_tJet_PosLep", &m_tJet_PosLep_phipair);
  m_reader->AddVariable("PhiPair_gJet_PosLep", &m_gJet_PosLep_phipair);
  m_reader->AddVariable("PhiPair_HcJet_NegLep", &m_HcJet_NegLep_phipair);
  m_reader->AddVariable("PhiPair_tHcJet_NegLep", &m_tHcJet_NegLep_phipair);
  m_reader->AddVariable("PhiPair_tJet_NegLep", &m_tJet_NegLep_phipair);
  m_reader->AddVariable("PhiPair_gJet_NegLep", &m_gJet_NegLep_phipair);
  m_reader->AddVariable("EPair_HcJet_PosLep", &m_HcJet_PosLep_epair);
  m_reader->AddVariable("EPair_tHcJet_PosLep", &m_tHcJet_PosLep_epair);
  m_reader->AddVariable("EPair_tJet_PosLep", &m_tJet_PosLep_epair);
  m_reader->AddVariable("EPair_gJet_PosLep", &m_gJet_PosLep_epair);
  m_reader->AddVariable("EPair_HcJet_NegLep", &m_HcJet_NegLep_epair);
  m_reader->AddVariable("EPair_tHcJet_NegLep", &m_tHcJet_NegLep_epair);
  m_reader->AddVariable("EPair_tJet_NegLep", &m_tJet_NegLep_epair);
  m_reader->AddVariable("EPair_gJet_NegLep", &m_gJet_NegLep_epair);


  m_reader->BookMVA( "BDT method", m_WeightsFile );

}

void BNWReconstruction::GetBDTWeight() {


  float weight =-9999999;

  ResetBDTVariables();

  for (int Hc = 0; Hc < 4; ++Hc) {
    for (int tHc = 0; tHc < 4; ++tHc) {
      if (tHc == Hc) continue;
      for (int t = 0; t < 4; ++t) {
        if (t == Hc || t == tHc) continue;
        for (int g = 0; g < 4; ++g) {
          if (g == Hc || g == tHc || g == t) continue;

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

          weight = m_reader->EvaluateMVA( "BDT method" );

          // std::cout << "Weight: "<< weight << " tHc " << tHc << " t " << t << std::endl;
          if ( weight > m_maxWeight ) {

            m_maxWeight = weight;

            m_HcJet_index  = Hc;
            m_tHcJet_index = tHc;
            m_tJet_index   = t;
            m_gJet_index   = g;

          }
          // std::cout << "Max Weight: "<< m_maxWeight << " m_tHcJet_index " << m_tHcJet_index << " m_tJet_index " << m_tJet_index << std::endl;


        }
      }
    }
  }


  return;


}

void BNWReconstruction::ResetBDTVariables() {
  m_maxWeight = -999999;
  m_HcJet_index  = -9;
  m_tHcJet_index = -9;
  m_tJet_index   = -9;
  m_gJet_index   = -9;

}

void BNWReconstruction::BookNewTree() {
    //---------------
  // Running NW on all combinations of b-partons (i.e. use v_bQuarksLorentz)
  //---------------

  m_outputTree = (TTree*) m_inputTree->CloneTree(0);

  m_outputTree->Branch("isReconstructed", &isReconstructed, "isReconstructed/I");

  m_outputTree->Branch("MaximumNeutrinoWeight", &maxNeutrinoWeight, "MaximumNeutrinoWeight/F");

  m_outputTree->Branch("MaxNW_top_Lorentz", "TLorentzVector", &maxNW_top);
  m_outputTree->Branch("MaxNW_top_pt", &maxNW_top_pt, "MaxNW_top_pt/F");
  m_outputTree->Branch("MaxNW_top_eta", &maxNW_top_eta, "MaxNW_top_eta/F");
  m_outputTree->Branch("MaxNW_top_phi", &maxNW_top_phi, "MaxNW_top_phi/F");
  m_outputTree->Branch("MaxNW_top_E", &maxNW_top_E, "MaxNW_top_E/F");
  m_outputTree->Branch("MaxNW_top_M", &maxNW_top_M, "MaxNW_top_M/F");

  m_outputTree->Branch("MaxNW_tbar_Lorentz", "TLorentzVector", &maxNW_tbar);
  m_outputTree->Branch("MaxNW_tbar_pt", &maxNW_tbar_pt, "MaxNW_tbar_pt/F");
  m_outputTree->Branch("MaxNW_tbar_eta", &maxNW_tbar_eta, "MaxNW_tbar_eta/F");
  m_outputTree->Branch("MaxNW_tbar_phi", &maxNW_tbar_phi, "MaxNW_tbar_phi/F");
  m_outputTree->Branch("MaxNW_tbar_E", &maxNW_tbar_E, "MaxNW_tbar_E/F");
  m_outputTree->Branch("MaxNW_tbar_M", &maxNW_tbar_M, "MaxNW_tbar_M/F");

  m_outputTree->Branch("MaxNW_b_Lorentz", "TLorentzVector", &maxNW_b);
  m_outputTree->Branch("MaxNW_b_Index", &maxNW_b_index, "MaxNW_b_Index/I");
  m_outputTree->Branch("MaxNW_b_pt", &maxNW_b_pt, "MaxNW_b_pt/F");
  m_outputTree->Branch("MaxNW_b_eta", &maxNW_b_eta, "MaxNW_b_eta/F");
  m_outputTree->Branch("MaxNW_b_phi", &maxNW_b_phi, "MaxNW_b_phi/F");
  m_outputTree->Branch("MaxNW_b_E", &maxNW_b_E, "MaxNW_b_E/F");
  m_outputTree->Branch("MaxNW_b_M", &maxNW_b_M, "MaxNW_b_M/F");

  m_outputTree->Branch("MaxNW_bbar_Lorentz", "TLorentzVector", &maxNW_bbar);
  m_outputTree->Branch("MaxNW_bbar_Index", &maxNW_bbar_index, "MaxNW_bbar_Index/I");
  m_outputTree->Branch("MaxNW_bbar_pt", &maxNW_bbar_pt, "MaxNW_bbar_pt/F");
  m_outputTree->Branch("MaxNW_bbar_eta", &maxNW_bbar_eta, "MaxNW_bbar_eta/F");
  m_outputTree->Branch("MaxNW_bbar_phi", &maxNW_bbar_phi, "MaxNW_bbar_phi/F");
  m_outputTree->Branch("MaxNW_bbar_E", &maxNW_bbar_E, "MaxNW_bbar_E/F");
  m_outputTree->Branch("MaxNW_bbar_M", &maxNW_bbar_M, "MaxNW_bbar_M/F");

  m_outputTree->Branch("MaxNW_nu_Lorentz", "TLorentzVector", &maxNW_nu);
  m_outputTree->Branch("MaxNW_nu_pt", &maxNW_nu_pt, "MaxNW_nu_pt/F");
  m_outputTree->Branch("MaxNW_nu_eta", &maxNW_nu_eta, "MaxNW_nu_eta/F");
  m_outputTree->Branch("MaxNW_nu_phi", &maxNW_nu_phi, "MaxNW_nu_phi/F");
  m_outputTree->Branch("MaxNW_nu_E", &maxNW_nu_E, "MaxNW_nu_E/F");
  m_outputTree->Branch("MaxNW_nu_M", &maxNW_nu_M, "MaxNW_nu_M/F");

  m_outputTree->Branch("MaxNW_nubar_Lorentz", "TLorentzVector", &maxNW_nubar);
  m_outputTree->Branch("MaxNW_nubar_pt", &maxNW_nubar_pt, "MaxNW_nubar_pt/F");
  m_outputTree->Branch("MaxNW_nubar_eta", &maxNW_nubar_eta, "MaxNW_nubar_eta/F");
  m_outputTree->Branch("MaxNW_nubar_phi", &maxNW_nubar_phi, "MaxNW_nubar_phi/F");
  m_outputTree->Branch("MaxNW_nubar_E", &maxNW_nubar_E, "MaxNW_nubar_E/F");
  m_outputTree->Branch("MaxNW_nubar_M", &maxNW_nubar_M, "MaxNW_nubar_M/F");

  m_outputTree->Branch("MaxNW_btH_Index", &maxNW_btH_index, "MaxNW_btH_Index/I");
  m_outputTree->Branch("MaxNW_bt_Index", &maxNW_bt_index, "MaxNW_bt_Index/I");

  m_outputTree->Branch("MaxNW_bH_Index", &maxNW_bH_index, "MaxNW_bH_Index/I");
  m_outputTree->Branch("MaxNW_bg_Index", &maxNW_bg_index, "MaxNW_bg_Index/I");
  m_outputTree->Branch("MaxBDTWeight"  , &MaxBDTWeight, "MaxBDTWeight/F");

}

void BNWReconstruction::Run() {

  long int allEntries = m_inputTree->GetEntries();

  GetBranches();

  SetTMVAReader();

  for ( int evt = 0; evt < allEntries; evt++ ) {

    if (evt%500 == 0) std::cout << "Event " << evt << " / " << allEntries << std::endl;
    InitialiseEventVariables();
    m_inputTree->GetEntry(evt);

    if (!passSelection) continue;

    GetBDTWeight(); // this assign the BDT reco variables

    double met_ex        = delphes_etmiss->Px();
    double met_ey        = delphes_etmiss->Py();
    double met_phi       = delphes_etmiss->Phi();

    NeutrinoWeighter nuW, nuW2;

    // std::cout << "PosLep Pt: " << delphes_poslep->Pt() << " NegLep Pt: " << delphes_neglep->Pt() << " Jet 1 Pt " << v_SelectedJets_Lorentz.at(m_tHcJet_index)->Pt() << "Jet 2 Pt" << v_SelectedJets_Lorentz.at(m_tJet_index)->Pt() << std::endl;

    double m_MaxWeight_1   = nuW.Reconstruct(*delphes_poslep, *delphes_neglep,  *v_SelectedJets_Lorentz.at(m_tHcJet_index),  *v_SelectedJets_Lorentz.at(m_tJet_index), met_ex, met_ey, met_phi);
    double m_MaxWeight_2   = nuW2.Reconstruct(*delphes_poslep, *delphes_neglep,  *v_SelectedJets_Lorentz.at(m_tJet_index),  *v_SelectedJets_Lorentz.at(m_tHcJet_index),  met_ex, met_ey, met_phi);

    // std::cout << "event: " << evt << " m_MaxWeight_1 " << m_MaxWeight_1 << " m_MaxWeight_2 " << m_MaxWeight_2 <<" Max BDT Weight: " << m_maxWeight << " m_tHcJet_index " << m_tHcJet_index << " m_tJet_index " << m_tJet_index << std::endl;

    maxNW_bH_index = m_HcJet_index;
    maxNW_bg_index = m_gJet_index;
    maxNW_btH_index =m_tHcJet_index;

    maxNW_bt_index  = m_tJet_index;
    MaxBDTWeight  = m_maxWeight;

    if (m_MaxWeight_1 > m_MaxWeight_2) {
      maxNeutrinoWeight = m_MaxWeight_1;
      maxNW_b_index     = m_tHcJet_index;
      maxNW_bbar_index  = m_tJet_index;

      maxNW_top        = nuW.GetTop();
      maxNW_tbar       = nuW.GetTbar();
      maxNW_b          = nuW.GetB();
      maxNW_bbar       = nuW.GetBbar();
      maxNW_nu         = nuW.GetNu();
      maxNW_nubar      = nuW.GetNubar();
    } else {
      maxNeutrinoWeight = m_MaxWeight_2;
      maxNW_b_index     = m_tJet_index;
      maxNW_bbar_index  = m_tHcJet_index;

      maxNW_top        = nuW2.GetTop();
      maxNW_tbar       = nuW2.GetTbar();
      maxNW_b          = nuW2.GetB();
      maxNW_bbar       = nuW2.GetBbar();
      maxNW_nu         = nuW2.GetNu();
      maxNW_nubar      = nuW2.GetNubar();
    }


    if (maxNeutrinoWeight > 0 ) isReconstructed = 1;
    else isReconstructed = 0;

    // std::cout << "Top Pt " << maxNW_top.Pt() << " Tbar PT " << maxNW_tbar.Pt() <<std::endl;
    maxNW_top_pt  = maxNW_top.Pt();
    maxNW_top_eta = maxNW_top.Eta();
    maxNW_top_phi = maxNW_top.Phi();
    maxNW_top_E   = maxNW_top.E();
    maxNW_top_M   = maxNW_top.M();

    maxNW_tbar_pt  = maxNW_tbar.Pt();
    maxNW_tbar_eta = maxNW_tbar.Eta();
    maxNW_tbar_phi = maxNW_tbar.Phi();
    maxNW_tbar_E   = maxNW_tbar.E();
    maxNW_tbar_M   = maxNW_tbar.M();

    maxNW_b_pt  = maxNW_b.Pt();
    maxNW_b_eta = maxNW_b.Eta();
    maxNW_b_phi = maxNW_b.Phi();
    maxNW_b_E   = maxNW_b.E();
    maxNW_b_M   = maxNW_b.M();

    maxNW_bbar_pt  = maxNW_bbar.Pt();
    maxNW_bbar_eta = maxNW_bbar.Eta();
    maxNW_bbar_phi = maxNW_bbar.Phi();
    maxNW_bbar_E   = maxNW_bbar.E();
    maxNW_bbar_M   = maxNW_bbar.M();

    maxNW_nu_pt  = maxNW_nu.Pt();
    maxNW_nu_eta = maxNW_nu.Eta();
    maxNW_nu_phi = maxNW_nu.Phi();
    maxNW_nu_E   = maxNW_nu.E();
    maxNW_nu_M   = maxNW_nu.M();

    maxNW_nubar_pt  = maxNW_nubar.Pt();
    maxNW_nubar_eta = maxNW_nubar.Eta();
    maxNW_nubar_phi = maxNW_nubar.Phi();
    maxNW_nubar_E   = maxNW_nubar.E();
    maxNW_nubar_M   = maxNW_nubar.M();

    m_outputTree->Fill();



  }

  m_outputTree->Write();
  m_outputFile->Close();

}




float BNWReconstruction::DPtPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->Pt() - b->Pt());
}

float BNWReconstruction::DEtaPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->Eta() - b->Eta());
}

float BNWReconstruction::DPhiPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->DeltaPhi(*b));
}

float BNWReconstruction::DMPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->M() - b->M());
}

float BNWReconstruction::DRPair(TLorentzVector* a, TLorentzVector* b) {
  return (a->DeltaR(*b));
}

float BNWReconstruction::PtPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).Pt();
}

float BNWReconstruction::MinvPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).M();
}

float BNWReconstruction::EtaPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).Eta();
}

float BNWReconstruction::PhiPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).Phi();
}

float BNWReconstruction::EPair(TLorentzVector* a, TLorentzVector* b) {
 return ((*a)+(*b)).E();
}
