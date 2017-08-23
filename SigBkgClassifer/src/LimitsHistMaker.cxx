#include "include/LimitsHistMaker.h"

LimitsHistMaker::LimitsHistMaker(std::string inputFileName, std::string outputFileName, std::string typeBDT, bool haveBDTWeightCut, float bdtWeightCutValue, float normalisation) :
m_inputFileName(inputFileName), m_outputFileName(outputFileName),
m_haveBDTWeightCut(haveBDTWeightCut), m_BDTWeightCutValue(bdtWeightCutValue),
m_norm(normalisation) {
  m_typeBDT = typeBDT;
  m_WeightsFile = "weights/TMVAClassification_BDT_"+typeBDT+".weights.xml";

  m_inputFile = new TFile(inputFileName.c_str(),"READ");
  m_inputTree = (TTree*)m_inputFile->Get("RecoResults");
  m_outputFile = new TFile(outputFileName.c_str(),"RECREATE");
}

void LimitsHistMaker::BookHists() {

  // add variables for sig/bkg discrimination BDT
  // h_numJets = new TH1F("NumberOfJets", "NumberOfJets", 9, 3, 12);
  h_maxRecoBDTWeight = new TH1F("maxRecoBDTWeight", "maxRecoBDTWeight", 16, -0.125, 0.275);
  // h_minv_bLeading_ClosestTop = new TH1F("Minv_bLeading_ClosestTop", "Minv_bLeading_ClosestTop", 120, 0., 1200.);
  h_minv_bH_ClosestTop = new TH1F("Minv_bH_ClosestTop", "Minv_bH_ClosestTop", 11, 175., 450.);
  // h_minv_minDRJets = new TH1F("Minv_minDRJets", "Minv_minDRJets", 240, 0., 1200.);

  int nbins_Ht  = 0;
  float xmin_Ht = 0;
  float xmax_Ht = 0;
  if (m_inputFileName.find("m300") != std::string::npos) {
    nbins_Ht = 15;
    xmin_Ht  = 200.;
    xmax_Ht  = 950.;
  }
  else if (m_inputFileName.find("m400") != std::string::npos) {
    nbins_Ht = 15;
    xmin_Ht  = 250.;
    xmax_Ht  = 1000.;
  }
  else if (m_inputFileName.find("m500") != std::string::npos) {
    nbins_Ht = 15;
    xmin_Ht  = 250.;
    xmax_Ht  = 1000.;
  }
  else if (m_inputFileName.find("m600") != std::string::npos) {
    nbins_Ht = 20;
    xmin_Ht  = 200.;
    xmax_Ht  = 1200.;
  }
  else if (m_inputFileName.find("m700") != std::string::npos) {
    nbins_Ht = 20;
    xmin_Ht  = 300.;
    xmax_Ht  = 1300.;
  }
  else if (m_inputFileName.find("m800") != std::string::npos) {
    nbins_Ht = 22;
    xmin_Ht  = 300.;
    xmax_Ht  = 1400.;
  }
  else if (m_inputFileName.find("m900") != std::string::npos) {
    nbins_Ht = 22;
    xmin_Ht  = 300.;
    xmax_Ht  = 1400.;
  }
  else if (m_inputFileName.find("m900") != std::string::npos) {
    nbins_Ht = 22;
    xmin_Ht  = 300.;
    xmax_Ht  = 1400.;
  }
  h_Ht = new TH1F("Ht", "Ht", nbins_Ht, xmin_Ht, xmax_Ht);

  h_minMinv_poslep_b = new TH1F("MinMinv_poslep_b", "MinMinv_poslep_b", 9, 20., 200.);
  h_minMinv_neglep_b = new TH1F("MinMinv_neglep_b", "MinMinv_neglep_b", 10, 20., 220.);
  // h_pt_leadingJet = new TH1F("Pt_LeadingJet", "Pt_LeadingJet", 100, 0., 1000.);
  h_pt_bH = new TH1F("Pt_bH", "Pt_bH", 15, 20., 320.);

  int nbins_minvHiggs  = 0;
  float xmin_minvHiggs = 0;
  float xmax_minvHiggs = 0;
  if (m_inputFileName.find("m300") != std::string::npos) {
    nbins_minvHiggs = 12;
    xmin_minvHiggs  = 200.;
    xmax_minvHiggs  = 500.;
  }
  else if (m_inputFileName.find("m400") != std::string::npos) {
    nbins_minvHiggs = 16;
    xmin_minvHiggs  = 200.;
    xmax_minvHiggs  = 600.;
  }
  else if (m_inputFileName.find("m500") != std::string::npos) {
    nbins_minvHiggs = 24;
    xmin_minvHiggs  = 200.;
    xmax_minvHiggs  = 800.;
  }
  else if (m_inputFileName.find("m600") != std::string::npos) {
    nbins_minvHiggs = 32;
    xmin_minvHiggs  = 200.;
    xmax_minvHiggs  = 1000.;
  }
  else if (m_inputFileName.find("m700") != std::string::npos) {
    nbins_minvHiggs = 40;
    xmin_minvHiggs  = 300.;
    xmax_minvHiggs  = 1300.;
  }
  else if (m_inputFileName.find("m800") != std::string::npos) {
    nbins_minvHiggs = 44;
    xmin_minvHiggs  = 300.;
    xmax_minvHiggs  = 1400.;
  }
  else if (m_inputFileName.find("m900") != std::string::npos) {
    nbins_minvHiggs = 44;
    xmin_minvHiggs  = 300.;
    xmax_minvHiggs  = 1400.;
  }
  h_minv_Higgs = new TH1F("Minv_Higgs", "Minv_Higgs", nbins_minvHiggs, xmin_minvHiggs, xmax_minvHiggs);

  h_dE_NonTopJets = new TH1F("DE_NonTopJets", "DE_NonTopJets", 13, -900., 400.);
  h_dR_Tops = new TH1F("DR_Tops", "DR_Tops", 12, 0., 6.);
  h_dR_bH_tH = new TH1F("DR_bH_tH", "DR_bH_tH", 9, 0., 4.5);
  h_dR_bLeading_TopLeading = new TH1F("DR_bLeading_TopLeading", "DR_bLeading_TopLeading", 10, 0., 5.);
  h_centrality = new TH1F("Centrality", "Centrality", 8, 0.2, 1.);
  h_dR_tH_lt = new TH1F("DR_tH_lt", "DR_tH_lt", 12, 0., 5.);
  h_dR_lep_maxDR_t_bH = new TH1F("DR_lep_maxDR_t_bH", "DR_lep_maxDR_t_bH", 11, 0., 5.5);
  h_costheta_lepton_jet = new TH1F("costheta_lepton_jet", "costheta_lepton_jet", 15, -1., 0.5);
  h_eta_bH = new TH1F("EtabH", "EtabH", 14, -3., 2.6);
  h_eta_btH = new TH1F("EtabtH", "EtabtH", 13, -2.6, 2.6);
  h_dR_bH_btH = new TH1F("DR_bH_btH", "DR_bH_btH", 8, 0., 4.);
  h_dR_bg_bt = new TH1F("DR_bg_bt", "DR_bg_bt", 10, 0., 5.);
  h_dEta_bg_bt = new TH1F("DeltaEta_bg_bt", "DeltaEta_bg_bt", 18, -4.5, 4.5);
  h_dPhi_btH_bH = new TH1F("DeltaPhi_btH_bH", "DeltaPhi_btH_bH", 14, -3.5, 3.5);
  h_dPhi_bt_bH = new TH1F("DeltaPhi_bt_bH", "DeltaPhi_bt_bH", 14, -3.5, 3.5);
  h_dPhi_bt_btH = new TH1F("DeltaPhi_bt_btH", "DeltaPhi_bt_btH", 14, -3.5, 3.5);
  h_DR_bg_bH = new TH1F("DR_bg_bH", "DR_bg_bH", 10, 0., 5.);
  h_dphi_poslep_neglep = new TH1F("DeltaPhi_poslep_neglep", "DeltaPhi_poslep_neglep", 14, -3.5, 3.5);
  h_dphi_ltH_bH = new TH1F("DeltaPhi_ltH_bH", "DeltaPhi_ltH_bH", 14, -3.5, 3.5);
  h_DR_ltH_bH = new TH1F("DR_ltH_bH", "DR_ltH_bH", 10, 0., 5.);
  h_dphi_ltH_bt = new TH1F("DeltaPhi_ltH_bt", "DeltaPhi_ltH_bt", 14, -3.5, 3.5);
  h_DR_ltH_bt = new TH1F("DR_ltH_bt", "DR_ltH_bt", 10, 0., 5.);
  return;
}

void LimitsHistMaker::WriteHists() {
  // h_numJets->Scale(m_norm);
  // h_numJets->Write();
  h_maxRecoBDTWeight->Scale(m_norm);
  h_maxRecoBDTWeight->Write();
  // h_minv_bLeading_ClosestTop->Scale(m_norm);
  // h_minv_bLeading_ClosestTop->Write();
  h_minv_bH_ClosestTop->Scale(m_norm);
  h_minv_bH_ClosestTop->Write();
  // h_minv_minDRJets->Scale(m_norm);
  // h_minv_minDRJets->Write();
  h_Ht->Scale(m_norm);
  h_Ht->Write();
  h_minMinv_poslep_b->Scale(m_norm);
  h_minMinv_poslep_b->Write();
  h_minMinv_neglep_b->Scale(m_norm);
  h_minMinv_neglep_b->Write();
  // h_pt_leadingJet->Scale(m_norm);
  // h_pt_leadingJet->Write();
  h_pt_bH->Scale(m_norm);
  h_pt_bH->Write();
  h_minv_Higgs->Scale(m_norm);
  h_minv_Higgs->Write();
  h_dE_NonTopJets->Scale(m_norm);
  h_dE_NonTopJets->Write();
  h_dR_Tops->Scale(m_norm);
  h_dR_Tops->Write();
  h_dR_bH_tH->Scale(m_norm);
  h_dR_bH_tH->Write();
  h_dR_bLeading_TopLeading->Scale(m_norm);
  h_dR_bLeading_TopLeading->Write();
  h_centrality->Scale(m_norm);
  h_centrality->Write();
  h_dR_tH_lt->Scale(m_norm);
  h_dR_tH_lt->Write();
  h_dR_lep_maxDR_t_bH->Scale(m_norm);
  h_dR_lep_maxDR_t_bH->Write();
  h_costheta_lepton_jet->Scale(m_norm);
  h_costheta_lepton_jet->Write();
  h_eta_bH->Scale(m_norm);
  h_eta_bH->Write();
  h_eta_btH->Scale(m_norm);
  h_eta_btH->Write();
  h_dR_bH_btH->Scale(m_norm);
  h_dR_bH_btH->Write();
  h_dR_bg_bt->Scale(m_norm);
  h_dR_bg_bt->Write();
  h_dEta_bg_bt->Scale(m_norm);
  h_dEta_bg_bt->Write();
  h_dPhi_btH_bH->Scale(m_norm);
  h_dPhi_btH_bH->Write();
  h_dPhi_bt_bH->Scale(m_norm);
  h_dPhi_bt_bH->Write();
  h_dPhi_bt_btH->Scale(m_norm);
  h_dPhi_bt_btH->Write();
  h_DR_bg_bH->Scale(m_norm);
  h_DR_bg_bH->Write();
  h_dphi_poslep_neglep->Scale(m_norm);
  h_dphi_poslep_neglep->Write();
  h_dphi_ltH_bH->Scale(m_norm);
  h_dphi_ltH_bH->Write();
  h_DR_ltH_bH->Scale(m_norm);
  h_DR_ltH_bH->Write();
  h_dphi_ltH_bt->Scale(m_norm);
  h_dphi_ltH_bt->Write();
  h_DR_ltH_bt->Scale(m_norm);
  h_DR_ltH_bt->Write();
  return;
}

void LimitsHistMaker::GetBranches() {

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

  m_inputTree->SetBranchAddress("isReconstructed", &isReconstructed);
  m_inputTree->SetBranchAddress("MaximumNeutrinoWeight", &maxNeutrinoWeight);

  maxNW_top = new TLorentzVector();
  m_inputTree->SetBranchAddress("MaxNW_top_Lorentz", &maxNW_top);
  m_inputTree->SetBranchAddress("MaxNW_top_pt", &maxNW_top_pt);
  m_inputTree->SetBranchAddress("MaxNW_top_eta", &maxNW_top_eta);
  m_inputTree->SetBranchAddress("MaxNW_top_phi", &maxNW_top_phi);
  m_inputTree->SetBranchAddress("MaxNW_top_E", &maxNW_top_E);
  m_inputTree->SetBranchAddress("MaxNW_top_M", &maxNW_top_M);

  maxNW_tbar = new TLorentzVector();
  m_inputTree->SetBranchAddress("MaxNW_tbar_Lorentz", &maxNW_tbar);
  m_inputTree->SetBranchAddress("MaxNW_tbar_pt", &maxNW_tbar_pt);
  m_inputTree->SetBranchAddress("MaxNW_tbar_eta", &maxNW_tbar_eta);
  m_inputTree->SetBranchAddress("MaxNW_tbar_phi", &maxNW_tbar_phi);
  m_inputTree->SetBranchAddress("MaxNW_tbar_E", &maxNW_tbar_E);
  m_inputTree->SetBranchAddress("MaxNW_tbar_M", &maxNW_tbar_M);

  maxNW_b = new TLorentzVector();
  m_inputTree->SetBranchAddress("MaxNW_b_Lorentz", &maxNW_b);
  m_inputTree->SetBranchAddress("MaxNW_b_Index", &maxNW_b_index);
  m_inputTree->SetBranchAddress("MaxNW_b_pt", &maxNW_b_pt);
  m_inputTree->SetBranchAddress("MaxNW_b_eta", &maxNW_b_eta);
  m_inputTree->SetBranchAddress("MaxNW_b_phi", &maxNW_b_phi);
  m_inputTree->SetBranchAddress("MaxNW_b_E", &maxNW_b_E);
  m_inputTree->SetBranchAddress("MaxNW_b_M", &maxNW_b_M);

  maxNW_bbar = new TLorentzVector();
  m_inputTree->SetBranchAddress("MaxNW_bbar_Lorentz", &maxNW_bbar);
  m_inputTree->SetBranchAddress("MaxNW_bbar_Index", &maxNW_bbar_index);
  m_inputTree->SetBranchAddress("MaxNW_bbar_pt", &maxNW_bbar_pt);
  m_inputTree->SetBranchAddress("MaxNW_bbar_eta", &maxNW_bbar_eta);
  m_inputTree->SetBranchAddress("MaxNW_bbar_phi", &maxNW_bbar_phi);
  m_inputTree->SetBranchAddress("MaxNW_bbar_E", &maxNW_bbar_E);
  m_inputTree->SetBranchAddress("MaxNW_bbar_M", &maxNW_bbar_M);

  maxNW_nu = new TLorentzVector();
  m_inputTree->SetBranchAddress("MaxNW_nu_Lorentz", &maxNW_nu);
  m_inputTree->SetBranchAddress("MaxNW_nu_pt", &maxNW_nu_pt);
  m_inputTree->SetBranchAddress("MaxNW_nu_eta", &maxNW_nu_eta);
  m_inputTree->SetBranchAddress("MaxNW_nu_phi", &maxNW_nu_phi);
  m_inputTree->SetBranchAddress("MaxNW_nu_E", &maxNW_nu_E);
  m_inputTree->SetBranchAddress("MaxNW_nu_M", &maxNW_nu_M);

  maxNW_nubar = new TLorentzVector();
  m_inputTree->SetBranchAddress("MaxNW_nubar_Lorentz", &maxNW_nubar);
  m_inputTree->SetBranchAddress("MaxNW_nubar_pt", &maxNW_nubar_pt);
  m_inputTree->SetBranchAddress("MaxNW_nubar_eta", &maxNW_nubar_eta);
  m_inputTree->SetBranchAddress("MaxNW_nubar_phi", &maxNW_nubar_phi);
  m_inputTree->SetBranchAddress("MaxNW_nubar_E", &maxNW_nubar_E);
  m_inputTree->SetBranchAddress("MaxNW_nubar_M", &maxNW_nubar_M);

  m_inputTree->SetBranchAddress("MaxNW_btH_Index", &maxNW_btH_index);
  m_inputTree->SetBranchAddress("MaxNW_bt_Index", &maxNW_bt_index);
  m_inputTree->SetBranchAddress("MaxNW_bH_Index", &maxNW_bH_index);
  m_inputTree->SetBranchAddress("MaxNW_bg_Index", &maxNW_bg_index);

  m_inputTree->SetBranchAddress("MaxBDTWeight", &maxRecoBDTWeight);

}

void LimitsHistMaker::SetTMVAReader() {

  m_reader = new TMVA::Reader("!Color:!Silent");

  // m_reader->AddSpectator("EventNumber");
  // m_reader->AddVariable("NumberOfJets", &m_njets);
  m_reader->AddVariable("MaxBDTWeight", &m_maxRecoBDTWeight);
  // m_reader->AddVariable("Minv_bLeading_ClosestTop", &m_minv_bleading_closesttop);
  m_reader->AddVariable("Minv_bH_ClosestTop", &m_minv_bH_closesttop);
  m_reader->AddVariable("Minv_minDRJets", &m_minv_minDRJets);
  //m_reader->AddVariable("Ht", &m_ht);
  m_reader->AddVariable("MinMinv_poslep_b", &m_minMinv_poslep_b);
  m_reader->AddVariable("MinMinv_neglep_b", &m_minMinv_neglep_b);
  // m_reader->AddVariable("Pt_LeadingJet", &m_leadingJetPt);
  m_reader->AddVariable("Pt_bH", &m_bHPt);
  m_reader->AddVariable("Minv_Higgs", &m_minv_Higgs);
  m_reader->AddVariable("DE_NonTopJets", &m_DE_nonTopJets);
  m_reader->AddVariable("DR_Tops", &m_DR_tops);
  m_reader->AddVariable("DR_bH_tH", &m_DR_bH_tH);
  m_reader->AddVariable("DR_bLeading_TopLeading", &m_DR_bleading_TopLeading);
  m_reader->AddVariable("Centrality", &m_centrality);
  m_reader->AddVariable("DR_tH_lt", &m_DR_tH_lt);
  m_reader->AddVariable("DR_lep_maxDR_t_bH", &m_DR_lep_maxDR_t_bH);
  m_reader->AddVariable("costheta_lepton_jet", &m_costheta_lepton_jet);
  m_reader->AddVariable("EtabH", &m_eta_bH);
  m_reader->AddVariable("EtabtH", &m_eta_btH);
  m_reader->AddVariable("DeltaEta_bg_bt", &m_deta_bg_bt);
  m_reader->AddVariable("DeltaPhi_btH_bH", &m_dphi_btH_bH);
  m_reader->AddVariable("DeltaPhi_bt_bH", &m_dphi_bt_bH);
  m_reader->AddVariable("DR_bH_btH", &m_DR_bH_btH);
  m_reader->AddVariable("DR_bg_bt", &m_DR_bg_bt);
  m_reader->AddVariable("DeltaPhi_bt_btH", &m_dphi_bt_btH);
  m_reader->AddVariable("DR_bg_bH", &m_DR_bg_bH);
  m_reader->AddVariable("DeltaPhi_poslep_neglep", &m_dphi_poslep_neglep);
  m_reader->AddVariable("DeltaPhi_ltH_bH", &m_dphi_ltH_bH);
  m_reader->AddVariable("DR_ltH_bH", &m_DR_ltH_bH);
  m_reader->AddVariable("DeltaPhi_ltH_bt", &m_dphi_ltH_bt);
  m_reader->AddVariable("DR_ltH_bt", &m_DR_ltH_bt);

  TString methodName = "BDT method";
  m_reader->BookMVA(methodName, m_WeightsFile);

  return;
}

void LimitsHistMaker::GetBDTResult() {

  m_BDTWeight = -9999999;

  topFromHiggs   = NULL;
  nonHiggsLepton = NULL;
  higgsLepton = NULL;

  if (maxNW_btH_index == maxNW_b_index) {
    topFromHiggs   = maxNW_top; // Hplus
    nonHiggsLepton = delphes_neglep;
    higgsLepton    = delphes_poslep;
  }
  else if (maxNW_btH_index == maxNW_bbar_index) {
    topFromHiggs   = maxNW_tbar; // Hmin
    nonHiggsLepton = delphes_poslep;
    higgsLepton    = delphes_neglep;
  }

  if (maxNW_btH_index == maxNW_b_index) m_reconstructedCharge = 1;
  else m_reconstructedCharge = -1;

  bFromTop         = v_SelectedJets_Lorentz.at(maxNW_bt_index);
  bFromHiggs       = v_SelectedJets_Lorentz.at(maxNW_bH_index);
  bFromTopAndHiggs = v_SelectedJets_Lorentz.at(maxNW_btH_index);
  bFromGluon       = v_SelectedJets_Lorentz.at(maxNW_bg_index);

  //---------------
  // variables for BDT
  //---------------
  m_maxRecoBDTWeight         = maxRecoBDTWeight;
  // m_njets                    = (float)numJets;
  // get invariant mass of b and top (bH or leading b)
  // m_minv_bleading_closesttop = GetMinv_b_ClosestTop(v_SelectedJets_Lorentz.at(0));
  m_minv_bH_closesttop       = GetMinv_b_ClosestTop(bFromHiggs);
  // invariant mass of closest jets
  m_minv_minDRJets           = GetMinv_JetsWithMinDR();
  // minimum invariant mass between lepton and b
  m_minMinv_poslep_b         = GetMinvMin_lep_b("positive");
  m_minMinv_neglep_b         = GetMinvMin_lep_b("negative");
  // pt or leading jet or bH jet
  // m_leadingJetPt             = v_SelectedJets_Lorentz.at(0)->Pt();
  m_bHPt                     = bFromHiggs->Pt();
  // b-tagged jet multiplicity
  m_nBTags                   = numBTags;
  // get invariant mass of bH and tH
  m_minv_Higgs               = GetMinvPair(topFromHiggs,bFromHiggs);
  // dE between non-top jets (not used in NW)
  m_DE_nonTopJets            = GetDE_NonTopJets();
  // dR between tops
  m_DR_tops                  = GetDR(maxNW_top, maxNW_tbar);
  // get DR between tH and bH
  m_DR_bH_tH                 = GetDR(v_SelectedJets_Lorentz.at(maxNW_bH_index), topFromHiggs);
  // centrality
  m_centrality               = GetCentrality();
  // dR between top from Higgs and non-Higgs lepton
  m_DR_tH_lt                 = GetDR(topFromHiggs,nonHiggsLepton);
  // costheta
  m_costheta_lepton_jet      = GetCosTheta();
  // eta, deltaEta and deltaPhi
  m_eta_bH                   = bFromHiggs->Eta();
  m_eta_btH                  = bFromTopAndHiggs->Eta();
  m_deta_bg_bt               = bFromGluon->Eta() - bFromTop->Eta();
  m_dphi_btH_bH              = bFromTopAndHiggs->DeltaPhi(*bFromHiggs);
  m_dphi_bt_bH               = bFromTop->DeltaPhi(*bFromHiggs);
  m_dphi_bt_btH              = bFromTop->DeltaPhi(*bFromTopAndHiggs);
  m_dphi_poslep_neglep       = delphes_poslep->DeltaPhi(*delphes_neglep);
  m_dphi_ltH_bH              = higgsLepton->DeltaPhi(*bFromHiggs);
  m_dphi_ltH_bt              = higgsLepton->DeltaPhi(*bFromTop);
  // deltaR
  m_DR_bH_btH                = GetDR(bFromHiggs, bFromTopAndHiggs);
  m_DR_bg_bt                 = GetDR(bFromGluon, bFromTop);
  m_DR_bg_bH                 = GetDR(bFromGluon, bFromHiggs);
  m_DR_ltH_bH                = GetDR(higgsLepton, bFromHiggs);
  m_DR_ltH_bt                = GetDR(higgsLepton, bFromTop);

  // get DR between leading top and leading b-jet
  TLorentzVector* bleading = v_SelectedJets_Lorentz.at(0);
  if (maxNW_top->Pt() > maxNW_tbar->Pt()) {
    m_DR_bleading_TopLeading = GetDR(v_SelectedJets_Lorentz.at(0), maxNW_top);
  }
  else {
    m_DR_bleading_TopLeading = GetDR(v_SelectedJets_Lorentz.at(0), maxNW_tbar);
  }

  // dR between the t-bH system with dRmax and the opposite sign lepton
  TLorentzVector systemHiggs;
  TLorentzVector* lep_opp = NULL;
  if (GetDR(bFromHiggs, maxNW_top) > GetDR(bFromHiggs, maxNW_tbar)) {
    systemHiggs = (*bFromHiggs) + (*maxNW_top);
    lep_opp     = delphes_neglep;
  }
  else {
    systemHiggs = (*bFromHiggs) + (*maxNW_tbar);
    lep_opp     = delphes_poslep;
  }
  if (lep_opp != NULL) m_DR_lep_maxDR_t_bH = lep_opp->DeltaR(systemHiggs);

  m_BDTWeight = m_reader->EvaluateMVA("BDT method");

  return;
}

void LimitsHistMaker::Run() {

  BookHists();
  GetBranches();
  SetTMVAReader();

  int n_pass = 0;
  int n_passBDT = 0;

  // run over events
  int nEntries = m_inputTree->GetEntries();
  for (int entry = 0; entry < nEntries; entry++) {
    m_inputTree->GetEntry(entry);

    m_maxRecoBDTWeight       = -999.;
    // m_njets = -999;
    // m_minv_bleading_closesttop = -999.;
    m_minv_bH_closesttop     = -999.;
    m_minv_minDRJets         = -999.;
    m_ht                     = -999.;
    m_minMinv_poslep_b       = -999.;
    m_minMinv_neglep_b       = -999.;
    m_bHPt                   = -999.;
    m_nBTags                 = -999;
    m_minv_Higgs             = -999.;
    m_DE_nonTopJets          = -999.;
    m_DR_tops                = -999.;
    m_DR_bH_tH               = -999.;
    m_centrality             = -999.;
    m_DR_tH_lt               = -999.;
    m_costheta_lepton_jet    = -999.;
    m_eta_bH                 = -999.;
    m_eta_btH                = -999.;
    m_deta_bg_bt             = -999.;
    m_dphi_btH_bH            = -999.;
    m_dphi_bt_bH             = -999.;
    m_dphi_bt_btH            = -999.;
    m_dphi_poslep_neglep     = -999.;
    m_dphi_ltH_bH            = -999.;
    m_dphi_ltH_bt            = -999.;
    m_DR_bH_btH              = -999.;
    m_DR_bg_bt               = -999.;
    m_DR_bleading_TopLeading = -999.;
    m_DR_lep_maxDR_t_bH      = -999.;
    m_DR_bg_bH               = -999.;
    m_DR_ltH_bH              = -999.;
    m_DR_ltH_bt              = -999.;

    if (isReconstructed == 1 && numBTags >= 4) {
      n_pass++;
      GetBDTResult(); // this assign the BDT reco variables
      bool passBDTWeightCut = false;
      if (m_haveBDTWeightCut) {
        passBDTWeightCut = (m_BDTWeight > m_BDTWeightCutValue);
      }
      else passBDTWeightCut = true; // case where there is no cut so all events will pass

      /*
      if (entry < 200) {
        std::cout << "Event " << entry << std::endl;
        std::cout << " Reconstructed charge = " << m_reconstructedCharge << std::endl;
        std::cout << " m_maxRecoBDTWeight = " << m_maxRecoBDTWeight << std::endl;
        std::cout << " m_minv_bH_closesttop = " << m_minv_bH_closesttop << std::endl;
        std::cout << " m_minv_minDRJets = " << m_minv_minDRJets << std::endl;
        std::cout << " m_minMinv_poslep_b = " << m_minMinv_poslep_b << std::endl;
        std::cout << " m_minMinv_neglep_b = " << m_minMinv_neglep_b << std::endl;
        std::cout << " m_bHPt = " << m_bHPt << std::endl;
        std::cout << " m_minv_Higgs = " << m_minv_Higgs << std::endl;
        std::cout << " m_DE_nonTopJets = " << m_DE_nonTopJets << std::endl;
        std::cout << " m_DR_tops = " << m_DR_tops << std::endl;
        std::cout << " m_DR_bH_tH = " << m_DR_bH_tH << std::endl;
        std::cout << " m_DR_bleading_TopLeading = " << m_DR_bleading_TopLeading << std::endl;
        std::cout << " m_centrality = " << m_centrality << std::endl;
        std::cout << " m_DR_tH_lt = " << m_DR_tH_lt << std::endl;
        std::cout << " m_DR_lep_maxDR_t_bH = " << m_DR_lep_maxDR_t_bH << std::endl;
        std::cout << " m_costheta_lepton_jet = " << m_costheta_lepton_jet << std::endl;
        std::cout << " m_eta_bH = " << m_eta_bH << std::endl;
        std::cout << " m_eta_btH = " << m_eta_btH << std::endl;
        std::cout << " m_deta_bg_bt = " << m_deta_bg_bt << std::endl;
        std::cout << " m_dphi_btH_bH = " << m_dphi_btH_bH << std::endl;
        std::cout << " m_dphi_bt_bH = " << m_dphi_bt_bH << std::endl;
        std::cout << " m_DR_bH_btH = " << m_DR_bH_btH << std::endl;
        std::cout << " m_DR_bg_bt = " << m_DR_bg_bt << std::endl;
        std::cout << " m_dphi_bt_btH = " << m_dphi_bt_btH << std::endl;
        std::cout << " !! m_BDTWeight = " << m_BDTWeight << std::endl;
        std::cout << " --------------" << std::endl;
      }
      */

      if (passBDTWeightCut) {
        n_passBDT++;

        // ht
        m_ht = GetHt();

        // h_numJets->Fill(numJets);
        h_maxRecoBDTWeight->Fill(m_maxRecoBDTWeight);
        // h_minv_bLeading_ClosestTop->Fill(m_minv_bleading_closesttop);
        h_minv_bH_ClosestTop->Fill(m_minv_bH_closesttop);
        //h_minv_minDRJets->Fill(m_minv_minDRJets);
        h_Ht->Fill(m_ht);
        h_minMinv_poslep_b->Fill(m_minMinv_poslep_b);
        h_minMinv_neglep_b->Fill(m_minMinv_neglep_b);
        // h_pt_leadingJet->Fill(m_leadingJetPt);
        h_pt_bH->Fill(m_bHPt);
        h_minv_Higgs->Fill(m_minv_Higgs);
        h_dE_NonTopJets->Fill(m_DE_nonTopJets);
        h_dR_Tops->Fill(m_DR_tops);
        h_dR_bH_tH->Fill(m_DR_bH_tH);
        h_dR_bLeading_TopLeading->Fill(m_DR_bleading_TopLeading);
        h_centrality->Fill(m_centrality);
        h_dR_tH_lt->Fill(m_DR_tH_lt);
        h_dR_lep_maxDR_t_bH->Fill(m_DR_lep_maxDR_t_bH);
        h_costheta_lepton_jet->Fill(m_costheta_lepton_jet);
        h_eta_bH->Fill(m_eta_bH);
        h_eta_btH->Fill(m_eta_btH);
        h_dR_bH_btH->Fill(m_DR_bH_btH);
        h_dR_bg_bt->Fill(m_DR_bg_bt);
        h_dEta_bg_bt->Fill(m_deta_bg_bt);
        h_dPhi_btH_bH->Fill(m_dphi_btH_bH);
        h_dPhi_bt_bH->Fill(m_dphi_bt_bH);
        h_dPhi_bt_btH->Fill(m_dphi_bt_btH);
        h_DR_bg_bH->Fill(m_DR_bg_bH);
        h_dphi_poslep_neglep->Fill(m_dphi_poslep_neglep);
        h_dphi_ltH_bH->Fill(m_dphi_ltH_bH);
        h_DR_ltH_bH->Fill(m_DR_ltH_bH);
        h_dphi_ltH_bt->Fill(m_dphi_ltH_bt);
        h_DR_ltH_bt->Fill(m_DR_ltH_bt);
      }
    }
  }

  std::cout << "Normalise histograms by " << m_norm << " and write to file. " << std::endl;
  std::cout << "Normalised total entries before cuts: " << nEntries*m_norm << " from " << nEntries << std::endl;
  std::cout << "Normalised total events after cuts:   " << n_passBDT*m_norm << " from " << n_passBDT << std::endl;

  WriteHists();
  DeleteHists();

  m_outputFile->Close();
  m_inputFile->Close();
  DeleteChainVariables();
  delete m_reader;
  nEntries = 0;
  n_passBDT = 0;
  m_norm = 1;

  return;
}


float LimitsHistMaker::GetMinv_b_ClosestTop(TLorentzVector* b) {

  if (GetDR(b, maxNW_top) < GetDR(b, maxNW_tbar)) {
    return GetMinvPair(maxNW_top, b);
  }
  else {
    return GetMinvPair(maxNW_tbar, b);
  }

  // should not get to here
  return -999999.;
}

float LimitsHistMaker::GetMinv_JetsWithMinDR() {

  float minDR = 999999.;
  float minv = -999.;

  for (auto jet1 : v_SelectedJets_Lorentz) {
    for (auto jet2 : v_SelectedJets_Lorentz) {
      if ((*jet1) == (*jet2)) continue;

      if (GetDR(jet1, jet2) < minDR) {
        minDR = GetDR(jet1, jet2);
        minv = GetMinvPair(jet1, jet2);
      }
    }
  }

  return minv;
}

float LimitsHistMaker::GetCosTheta() {

  TLorentzVector higgs_system, *lep, *jet;
  TVector3 lepvec, jetvec;

  if (m_reconstructedCharge == 1 ){
    higgs_system = *v_SelectedJets_Lorentz.at(maxNW_bH_index) + *maxNW_top;
    lep = delphes_poslep;
  }
  else {
    higgs_system = *v_SelectedJets_Lorentz.at(maxNW_bH_index) + *maxNW_tbar;
    lep = delphes_neglep;
  }

  jet = v_SelectedJets_Lorentz.at(maxNW_bH_index);

  lep->Boost(- higgs_system.BoostVector());
  jet->Boost(- higgs_system.BoostVector());
  lepvec = lep->Vect();
  jetvec = jet->Vect();

  return (lepvec.Dot(jetvec))/(lepvec.Mag()*jetvec.Mag());

}

float LimitsHistMaker::GetHt() {

  float ht = 0;

  for (auto jet : v_SelectedJets_Lorentz) ht += jet->Pt();
  ht += delphes_poslep->Pt();
  ht += delphes_neglep->Pt();

  return ht;
}

float LimitsHistMaker::GetMinvMin_lep_b(std::string charge) {

  float minv = 999999.;

  for (auto jet : v_SelectedJets_Lorentz) {
    if (minv > GetMinvPair(jet,delphes_poslep) && charge.compare("positive") == 0) {
      minv = GetMinvPair(jet,delphes_poslep);
    }
    else if (minv > GetMinvPair(jet,delphes_neglep) && charge.compare("negative") == 0) {
      minv = GetMinvPair(jet,delphes_neglep);
    }
  }

  return minv;
}

float LimitsHistMaker::GetDE_NonTopJets() {
  TLorentzVector bH = *(v_SelectedJets_Lorentz.at(maxNW_bH_index));
  TLorentzVector bg = *(v_SelectedJets_Lorentz.at(maxNW_bg_index));
  return (bH-bg).E();
}

float LimitsHistMaker::GetCentrality() {

  float ht = GetHt();
  float e = 0.;
  for (auto jet : v_SelectedJets_Lorentz) e += jet->E();
  e += delphes_poslep->E();
  e += delphes_neglep->E();

  return ht/e;
}

float LimitsHistMaker::GetDR(TLorentzVector* a, TLorentzVector* b) {
  return a->DeltaR(*b);
}

float LimitsHistMaker::GetMinvPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).M();
}

void LimitsHistMaker::DeleteChainVariables() {

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

  delete maxNW_top;
  delete maxNW_tbar;
  delete maxNW_b;
  delete maxNW_bbar;
  delete maxNW_nu;
  delete maxNW_nubar;

  return;
}

void LimitsHistMaker::DeleteHists() {

  // delete h_numJets;
  delete h_maxRecoBDTWeight;
  // delete h_minv_bLeading_ClosestTop;
  delete h_minv_bH_ClosestTop;
  // delete h_minv_minDRJets;
  delete h_Ht;
  delete h_minMinv_poslep_b;
  delete h_minMinv_neglep_b;
  // delete h_pt_leadingJet;
  delete h_pt_bH;
  delete h_minv_Higgs;
  delete h_dE_NonTopJets;
  delete h_dR_Tops;
  delete h_dR_bH_tH;
  delete h_dR_bLeading_TopLeading;
  delete h_centrality;
  delete h_dR_tH_lt;
  delete h_dR_lep_maxDR_t_bH;
  delete h_costheta_lepton_jet;
  delete h_eta_bH;
  delete h_eta_btH;
  delete h_dR_bH_btH;
  delete h_dR_bg_bt;
  delete h_dEta_bg_bt;
  delete h_dPhi_btH_bH;
  delete h_dPhi_bt_bH;
  delete h_dPhi_bt_btH;
}
