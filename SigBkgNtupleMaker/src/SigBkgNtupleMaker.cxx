#include "include/SigBkgNtupleMaker.h"

SigBkgNtupleMaker::SigBkgNtupleMaker(std::string inputFileName, std::string outputFileName) :
m_inputFileName(inputFileName), m_outputFileName(outputFileName) {

  m_inputFile = new TFile(inputFileName.c_str());
  m_inputTree = (TTree*)m_inputFile->Get("RecoResults");
  m_outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  BookNewTree();
}

void SigBkgNtupleMaker::BookNewTree() {

  m_outputTree = (TTree*)m_inputTree->CloneTree(0);

  // add variables for sig/bkg discrimination BDT
  m_outputTree->Branch("Minv_bLeading_ClosestTop", &m_minv_bleading_closesttop, "Minv_bLeading_ClosestTop/F");
  m_outputTree->Branch("Minv_bH_ClosestTop", &m_minv_bH_closesttop, "Minv_bH_ClosestTop/F");
  m_outputTree->Branch("Minv_minDRJets", &m_minv_minDRJets, "Minv_minDRJets/F");
  m_outputTree->Branch("Ht", &m_ht, "Ht/F");
  m_outputTree->Branch("MinMinv_poslep_b", &m_minMinv_poslep_b, "MinMinv_poslep_b/F");
  m_outputTree->Branch("MinMinv_neglep_b", &m_minMinv_neglep_b, "MinMinv_neglep_b/F");
  m_outputTree->Branch("BTaggedJetMultiplicity", &m_nBTags, "BTaggedJetMultiplicity/I");
  m_outputTree->Branch("Pt_LeadingJet", &m_leadingJetPt, "LeadingJetPt/F");
  m_outputTree->Branch("Pt_bH", &m_bHPt,"Pt_bH/F");
  m_outputTree->Branch("Minv_Higgs", &m_minv_Higgs, "Minv_Higgs/F");
  m_outputTree->Branch("DE_NonTopJets", &m_DE_nonTopJets, "DE_NonTopJets/F");
  m_outputTree->Branch("DR_Tops", &m_DR_tops, "DR_Tops/F");
  m_outputTree->Branch("DR_bH_tH", &m_DR_bH_tH, "DR_bH_tH/F");
  m_outputTree->Branch("DR_bLeading_TopLeading", &m_DR_bleading_TopLeading, "DR_bLeading_TopLeading/F");
  m_outputTree->Branch("Centrality", &m_centrality, "Centrality/F");
  m_outputTree->Branch("DR_tH_lt", &m_DR_tH_lt, "DR_tH_lt/F");
  m_outputTree->Branch("DR_lep_maxDR_t_bH", &m_DR_lep_maxDR_t_bH, "DR_lep_maxDR_t_bH/F");
  m_outputTree->Branch("costheta_lepton_jet", &m_costheta_lepton_jet, "costheta_lepton_jet/F");
  m_outputTree->Branch("EtabH", &m_eta_bH, "eta_bH/F");
  m_outputTree->Branch("EtabtH", &m_eta_btH, "eta_btH/F");

  m_outputTree->Branch("DR_bH_btH", &m_DR_bH_btH, "DR_bH_btH/F");
  m_outputTree->Branch("DR_bg_bt", &m_DR_bg_bt, "DR_bg_bt/F");
  m_outputTree->Branch("DR_bg_bH", &m_DR_bg_bH, "DR_bg_bH/F");

  m_outputTree->Branch("DeltaEta_bg_bt", &m_deta_bg_bt, "DeltaEta_bg_bt/F");
  m_outputTree->Branch("DeltaPhi_btH_bH", &m_dphi_btH_bH, "DeltaPhi_btH_bH/F");
  m_outputTree->Branch("DeltaPhi_bt_bH", &m_dphi_bt_bH, "DeltaPhi_bt_bH/F");
  m_outputTree->Branch("DeltaPhi_bt_btH", &m_dphi_bt_btH, "DeltaPhi_bt_btH/F");
  m_outputTree->Branch("DeltaPhi_poslep_neglep", &m_dphi_poslep_neglep, "DeltaPhi_poslep_neglep/F");

  m_outputTree->Branch("DeltaPhi_ltH_bH", &m_dphi_ltH_bH, "DeltaPhi_ltH_bH/F");
  m_outputTree->Branch("DR_ltH_bH", &m_DR_ltH_bH, "DR_ltH_bH/F");
  m_outputTree->Branch("DeltaPhi_ltH_bt", &m_dphi_ltH_bt, "DeltaPhi_ltH_bt/F");
  m_outputTree->Branch("DR_ltH_bt", &m_DR_ltH_bt, "DR_ltH_bt/F");

  return;
}

void SigBkgNtupleMaker::InitialiseEventVariables() {

  m_minv_bleading_closesttop = -999.;
  m_minv_bH_closesttop = -999.;
  m_minv_minDRJets = -999.;
  m_ht = -999.;
  m_minMinv_poslep_b = -999.;
  m_minMinv_neglep_b = -999.;
  m_nBTags = -99;
  m_leadingJetPt = -999.;
  m_bHPt = -999.;
  m_minv_Higgs = -999.;
  m_DE_nonTopJets = -999.;
  m_DR_tops = -999.;
  m_DR_bH_tH = -999.;
  m_DR_bleading_TopLeading = -999.;
  m_centrality = -999.;
  m_DR_tH_lt = -999.;
  m_DR_lep_maxDR_t_bH = -999.;
  m_costheta_lepton_jet = -999.;
  m_reconstructedCharge = 0;

  return;
}

void SigBkgNtupleMaker::GetBranches() {

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

  m_inputTree->SetBranchAddress("MaxBDTWeight", &maxBDTWeight);

}

void SigBkgNtupleMaker::Run() {

  GetBranches();

  // run over events
  int nEntries = m_inputTree->GetEntries();
  for (int entry = 0; entry < nEntries; entry++) {
    m_inputTree->GetEntry(entry);

    InitialiseEventVariables();

    if (isReconstructed && numBTags >= 2) {

      topFromHiggs = NULL;
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
      // get invariant mass of b and top (bH or leading b)
      m_minv_bleading_closesttop = GetMinv_b_ClosestTop(v_SelectedJets_Lorentz.at(0));
      m_minv_bH_closesttop       = GetMinv_b_ClosestTop(bFromHiggs);
      // invariant mass of closest jets
      m_minv_minDRJets           = GetMinv_JetsWithMinDR();
      // ht
      m_ht                       = GetHt();
      // minimum invariant mass between lepton and b
      m_minMinv_poslep_b         = GetMinvMin_lep_b("positive");
      m_minMinv_neglep_b         = GetMinvMin_lep_b("negative");
      // pt or leading jet or bH jet
      m_leadingJetPt             = v_SelectedJets_Lorentz.at(0)->Pt();
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
      // dR between top from Higgs (tH) and non-Higgs lepton (lt)
      m_DR_tH_lt                 = GetDR(topFromHiggs,nonHiggsLepton);
      // dPhi and dR between Higgs lepton (ltH) and bH
      m_dphi_ltH_bH              = higgsLepton->DeltaPhi(*bFromHiggs);
      m_DR_ltH_bH                = GetDR(higgsLepton,bFromHiggs);
      // dPhi and dR between Higgs lepton (ltH) and bt
      m_dphi_ltH_bt              = higgsLepton->DeltaPhi(*bFromTop);
      m_DR_ltH_bt                = GetDR(higgsLepton,bFromTop);

      m_costheta_lepton_jet      = GetCosTheta();

      m_eta_bH                   = bFromHiggs->Eta();
      m_eta_btH                  = bFromTopAndHiggs->Eta();
      m_deta_bg_bt               = bFromGluon->Eta() - bFromTop->Eta();
      m_dphi_btH_bH              = bFromTopAndHiggs->DeltaPhi(*bFromHiggs);
      m_dphi_bt_bH               = bFromTop->DeltaPhi(*bFromHiggs);
      m_dphi_poslep_neglep       = delphes_poslep->DeltaPhi(*delphes_neglep);

      m_DR_bH_btH                = GetDR(bFromHiggs, bFromTopAndHiggs);
      m_DR_bg_bt                 = GetDR(bFromGluon, bFromTop);
      m_DR_bg_bH                 = GetDR(bFromGluon, bFromHiggs);
      m_dphi_bt_btH              = bFromTop->DeltaPhi(*bFromTopAndHiggs);

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

      m_outputTree->Fill();
    }
  }

  m_outputTree->Write();

  m_outputFile->Close();
  m_inputFile->Close();
  DeleteChainVariables();

  return;
}

float SigBkgNtupleMaker::GetMinv_b_ClosestTop(TLorentzVector* b) {

  if (GetDR(b, maxNW_top) < GetDR(b, maxNW_tbar)) {
    return GetMinvPair(maxNW_top, b);
  }
  else {
    return GetMinvPair(maxNW_tbar, b);
  }

  // should not get to here
  return -999999.;
}

float SigBkgNtupleMaker::GetMinv_JetsWithMinDR() {

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

float SigBkgNtupleMaker::GetCosTheta() {

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

float SigBkgNtupleMaker::GetHt() {

  float ht = 0;

  for (auto jet : v_SelectedJets_Lorentz) ht += jet->Pt();
  ht += delphes_poslep->Pt();
  ht += delphes_neglep->Pt();

  return ht;
}

float SigBkgNtupleMaker::GetMinvMin_lep_b(std::string charge) {

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

float SigBkgNtupleMaker::GetDE_NonTopJets() {
  TLorentzVector bH = *(v_SelectedJets_Lorentz.at(maxNW_bH_index));
  TLorentzVector bg = *(v_SelectedJets_Lorentz.at(maxNW_bg_index));
  return (bH-bg).E();
}

float SigBkgNtupleMaker::GetCentrality() {

  float ht = GetHt();
  float e = 0.;
  for (auto jet : v_SelectedJets_Lorentz) e += jet->E();
  e += delphes_poslep->E();
  e += delphes_neglep->E();

  return ht/e;
}

float SigBkgNtupleMaker::GetDR(TLorentzVector* a, TLorentzVector* b) {
  return a->DeltaR(*b);
}

float SigBkgNtupleMaker::GetMinvPair(TLorentzVector* a, TLorentzVector* b) {
  return ((*a)+(*b)).M();
}

void SigBkgNtupleMaker::DeleteChainVariables() {

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
