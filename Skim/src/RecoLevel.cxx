#include "include/RecoLevel.h"

//---------------------
// PUBLIC VIRTUAL Set tree branches reading for filling
//---------------------
void RecoLevel::BookTreeBranches() {
  std::cout << "Creating RecoResults tree..." << std::endl;

  tree = new TTree("RecoResults","RecoResults");
  tree->Branch("EventNumber", &eventNum, "EventNumber/I");
  tree->Branch("PassSelection", &passSelection, "PassSelection/I");
  tree->Branch("NumberOfJets", &numJets, "NumberOfJets/I");
  tree->Branch("NumberBTaggedJets", &numBTags, "NumberBTaggedJets/I");
  tree->Branch("NumberMatchedJets", &numMatched, "NumberMatchedJets/I");

  tree->Branch("HasMatchSIG_bFromTopAndHiggs_bFromTop_bFromHiggs", &hasRequiredSIGMatches, "HasMatch_bFromTopAndHiggs_bFromTop_bFromHiggs/I");
  tree->Branch("HasMatchBKG_bFromTop1_bFromTop2", &hasRequiredBKGMatches, "HasMatch_bFromTop1_bFromTop2/I");

  //-------------
  // Truth variables
  //-------------
  if (m_sig) {
    tree->Branch("Truth_Hc_Lorentz", "TLorentzVector", &truth_Hc);
    tree->Branch("Truth_Hc_Charge", &truth_Hc_charge, "Truth_Hc_charge/I");
    tree->Branch("Truth_Hc_pt", &truth_Hc_pt, "Truth_Hc_pt/F");
    tree->Branch("Truth_Hc_eta",  &truth_Hc_eta, "Truth_Hc_eta/F");
    tree->Branch("Truth_Hc_phi", &truth_Hc_phi, "Truth_Hc_phi/F");
    tree->Branch("Truth_Hc_E", &truth_Hc_E, "Truth_Hc_E/F");
    tree->Branch("Truth_Hc_M", &truth_Hc_M, "Truth_Hc_M/F");
  }

  tree->Branch("Truth_top_Lorentz", "TLorentzVector", &truth_top);
  tree->Branch("Truth_top_pt", &truth_top_pt, "Truth_top_pt/F");
  tree->Branch("Truth_top_eta",  &truth_top_eta, "Truth_top_eta/F");
  tree->Branch("Truth_top_phi", &truth_top_phi, "Truth_top_phi/F");
  tree->Branch("Truth_top_E", &truth_top_E, "Truth_top_E/F");
  tree->Branch("Truth_top_M", &truth_top_M, "Truth_top_M/F");

  tree->Branch("Truth_tbar_Lorentz", "TLorentzVector", &truth_tbar);
  tree->Branch("Truth_tbar_pt", &truth_tbar_pt, "Truth_tbar_pt/F");
  tree->Branch("Truth_tbar_eta",  &truth_tbar_eta, "Truth_tbar_eta/F");
  tree->Branch("Truth_tbar_phi", &truth_tbar_phi, "Truth_tbar_phi/F");
  tree->Branch("Truth_tbar_E", &truth_tbar_E, "Truth_tbar_E/F");
  tree->Branch("Truth_tbar_M", &truth_tbar_M, "Truth_tbar_M/F");

  tree->Branch("Truth_PosLep_Lorentz", "TLorentzVector", &truth_poslep);
  tree->Branch("Truth_PosLep_pt", &truth_poslep_pt, "Truth_PosLep_pt/F");
  tree->Branch("Truth_PosLep_eta", &truth_poslep_eta, "Truth_PosLep_eta/F");
  tree->Branch("Truth_PosLep_phi", &truth_poslep_phi, "Truth_PosLep_phi/F");
  tree->Branch("Truth_PosLep_E", &truth_poslep_E, "Truth_PosLep_E/F");
  tree->Branch("Truth_PosLep_M", &truth_poslep_M, "Truth_PosLep_M/F");

  tree->Branch("Truth_NegLep_Lorentz", "TLorentzVector", &truth_neglep);
  tree->Branch("Truth_NegLep_pt", &truth_neglep_pt, "Truth_NegLep_pt/F");
  tree->Branch("Truth_NegLep_eta", &truth_neglep_eta, "Truth_NegLep_eta/F");
  tree->Branch("Truth_NegLep_phi", &truth_neglep_phi, "Truth_NegLep_phi/F");
  tree->Branch("Truth_NegLep_E", &truth_neglep_E, "Truth_NegLep_E/F");
  tree->Branch("Truth_NegLep_M", &truth_neglep_M, "Truth_NegLep_M/F");

  tree->Branch("Truth_nu_Lorentz", "TLorentzVector", &truth_nu);
  tree->Branch("Truth_nu_pt", &truth_nu_pt, "Truth_nu_pt/F");
  tree->Branch("Truth_nu_eta",  &truth_nu_eta, "Truth_nu_eta/F");
  tree->Branch("Truth_nu_phi", &truth_nu_phi, "Truth_nu_phi/F");
  tree->Branch("Truth_nu_E", &truth_nu_E, "Truth_nu_E/F");
  tree->Branch("Truth_nu_M", &truth_nu_M, "Truth_nu_M/F");

  tree->Branch("Truth_nubar_Lorentz", "TLorentzVector", &truth_nubar);
  tree->Branch("Truth_nubar_pt", &truth_nubar_pt, "Truth_nubar_pt/F");
  tree->Branch("Truth_nubar_eta",  &truth_nubar_eta, "Truth_nubar_eta/F");
  tree->Branch("Truth_nubar_phi", &truth_nubar_phi, "Truth_nubar_phi/F");
  tree->Branch("Truth_nubar_E", &truth_nubar_E, "Truth_nubar_E/F");
  tree->Branch("Truth_nubar_M", &truth_nubar_M, "Truth_nubar_M/F");

  tree->Branch("Truth_b0_Lorentz", "TLorentzVector", &v_bPartons_Lorentz.at(0));
  tree->Branch("Truth_b1_Lorentz", "TLorentzVector", &v_bPartons_Lorentz.at(1));
  tree->Branch("Truth_b2_Lorentz", "TLorentzVector", &v_bPartons_Lorentz.at(2));
  tree->Branch("Truth_b3_Lorentz", "TLorentzVector", &v_bPartons_Lorentz.at(3));

  tree->Branch("Truth_bPartons_pt", &v_bPartons_pt);
  tree->Branch("Truth_bPartons_eta", &v_bPartons_eta);
  tree->Branch("Truth_bPartons_phi", &v_bPartons_phi);
  tree->Branch("Truth_bPartons_E", &v_bPartons_E);
  tree->Branch("Truth_bPartons_M", &v_bPartons_M);
  tree->Branch("Truth_bPartons_Origin", &v_bPartons_Origin);

  //-------------
  // Reco information
  //-------------

  tree->Branch("Delphes_AllJets_pt", &v_AllJets_pt);
  tree->Branch("Delphes_AllJets_eta", &v_AllJets_eta);
  tree->Branch("Delphes_AllJets_phi", &v_AllJets_phi);
  tree->Branch("Delphes_AllJets_E", &v_AllJets_E);
  tree->Branch("Delphes_AllJets_M", &v_AllJets_M);
  tree->Branch("Delphes_AllJets_BTag", &v_AllJets_bTag);

  tree->Branch("Delphes_SelectedJet0_Lorentz", "TLorentzVector", &v_SelectedJets_Lorentz.at(0));
  tree->Branch("Delphes_SelectedJet1_Lorentz", "TLorentzVector", &v_SelectedJets_Lorentz.at(1));
  tree->Branch("Delphes_SelectedJet2_Lorentz", "TLorentzVector", &v_SelectedJets_Lorentz.at(2));
  tree->Branch("Delphes_SelectedJet3_Lorentz", "TLorentzVector", &v_SelectedJets_Lorentz.at(3));

  tree->Branch("Delphes_SelectedJets_pt", &v_SelectedJets_pt);
  tree->Branch("Delphes_SelectedJets_eta", &v_SelectedJets_eta);
  tree->Branch("Delphes_SelectedJets_phi", &v_SelectedJets_phi);
  tree->Branch("Delphes_SelectedJets_E", &v_SelectedJets_E);
  tree->Branch("Delphes_SelectedJets_M", &v_SelectedJets_M);
  tree->Branch("Delphes_SelectedJets_BTag", &v_SelectedJets_bTag);

  tree->Branch("MatchingInfo_IndexOfQuarkMatchedToJet", &v_bIndexMatchedToJet);
  tree->Branch("MatchingInfo_OriginOfQuarkMatchedToJet", &v_JetMatchedOrigin);

  tree->Branch("Delphes_PosLep_Lorentz", "TLorentzVector", &delphes_poslep);
  tree->Branch("Delphes_PosLep_pt", &delphes_poslep_pt, "Delphes_PosLep_pt/F");
  tree->Branch("Delphes_PosLep_eta", &delphes_poslep_eta, "Delphes_PosLep_eta/F");
  tree->Branch("Delphes_PosLep_phi", &delphes_poslep_phi, "Delphes_PosLep_phi/F");
  tree->Branch("Delphes_PosLep_E", &delphes_poslep_E, "Delphes_PosLep_E/F");
  tree->Branch("Delphes_PosLep_M", &delphes_poslep_M, "Delphes_PosLep_M/F");

  tree->Branch("Delphes_NegLep_Lorentz", "TLorentzVector", &delphes_neglep);
  tree->Branch("Delphes_NegLep_pt", &delphes_neglep_pt, "Delphes_NegLep_pt/F");
  tree->Branch("Delphes_NegLep_eta", &delphes_neglep_eta, "Delphes_NegLep_eta/F");
  tree->Branch("Delphes_NegLep_phi", &delphes_neglep_phi, "Delphes_NegLep_phi/F");
  tree->Branch("Delphes_NegLep_E", &delphes_neglep_E, "Delphes_NegLep_E/F");
  tree->Branch("Delphes_NegLep_M", &delphes_neglep_M, "Delphes_NegLep_M/F");

  tree->Branch("Delphes_ETmiss_Lorentz", "TLorentzVector", &delphes_etmiss);
  tree->Branch("Delphes_ETmiss_MET", &delphes_etmiss_met, "Delphes_ETmiss_MET/F");
  tree->Branch("Delphes_ETmiss_eta", &delphes_etmiss_eta, "Delphes_ETmiss_eta/F");
  tree->Branch("Delphes_ETmiss_phi", &delphes_etmiss_phi, "Delphes_ETmiss_phi/F");

  return;
}

//---------------------
// PRIVATE VIRTUAL Set tree variables to default value
//---------------------
void RecoLevel::InitialiseEventVariables() {

  numMatched         = -1; // if -1, means did not pass cuts
  passSelection      = 0; // 1 if passes cuts
  hasRequiredSIGMatches = 0; // check whether match to b-quarks from tops and from Higgs exist
  hasRequiredBKGMatches = 0; // check whether match to b-quarks from tops (no Higgs check for bkg)

  truth_Hc        = TLorentzVector();
  truth_Hc_charge = -99;
  truth_Hc_pt     = -99.;
  truth_Hc_eta    = -99.;
  truth_Hc_phi    = -99.;
  truth_Hc_E      = -99.;
  truth_Hc_M      = -99.;

  truth_top     = TLorentzVector();
  truth_top_pt  = -99.;
  truth_top_eta = -99.;
  truth_top_phi = -99.;
  truth_top_E   = -99.;
  truth_top_M   = -99.;

  truth_tbar     = TLorentzVector();
  truth_tbar_pt  = -99.;
  truth_tbar_eta = -99.;
  truth_tbar_phi = -99.;
  truth_tbar_E   = -99.;
  truth_tbar_M   = -99.;

  truth_poslep     = TLorentzVector();
  truth_poslep_pt  = -99.;
  truth_poslep_eta = -99.;
  truth_poslep_phi = -99.;
  truth_poslep_E   = -99.;
  truth_poslep_M   = -99.;

  truth_neglep     = TLorentzVector();
  truth_neglep_pt  = -99.;
  truth_neglep_eta = -99.;
  truth_neglep_phi = -99.;
  truth_neglep_E   = -99.;
  truth_neglep_M   = -99.;

  truth_nu     = TLorentzVector();
  truth_nu_pt  = -99.;
  truth_nu_eta = -99.;
  truth_nu_phi = -99.;
  truth_nu_E   = -99.;
  truth_nu_M   = -99.;

  truth_nubar     = TLorentzVector();
  truth_nubar_pt  = -99.;
  truth_nubar_eta = -99.;
  truth_nubar_phi = -99.;
  truth_nubar_E   = -99.;
  truth_nubar_M   = -99.;

  for (int i = 0; i < v_AllJetsLorentz.size(); i++) {
    v_AllJets_pt.push_back(-30.);
    v_AllJets_eta.push_back(-30.);
    v_AllJets_phi.push_back(-30.);
    v_AllJets_E.push_back(-30.);
    v_AllJets_M.push_back(-30.);
    v_AllJets_bTag.push_back(-1);
  }

  delphes_poslep     = TLorentzVector();
  delphes_poslep_pt  = -99.;
  delphes_poslep_eta = -99.;
  delphes_poslep_phi = -99.;
  delphes_poslep_E   = -99.;
  delphes_poslep_M   = -99.;

  delphes_neglep     = TLorentzVector();
  delphes_neglep_pt  = -99.;
  delphes_neglep_eta = -99.;
  delphes_neglep_phi = -99.;
  delphes_neglep_E   = -99.;
  delphes_neglep_M   = -99.;

  delphes_etmiss      = TLorentzVector();
  delphes_etmiss_met  = -99.;
  delphes_etmiss_eta  = -99.;
  delphes_etmiss_phi  = -99.;

  for (int i = 0; i < 4; i++) {
    v_bPartons_Lorentz.at(i) = TLorentzVector();
    v_bPartons_pt.at(i)      = -99.;
    v_bPartons_eta.at(i)     = -99.;
    v_bPartons_phi.at(i)     = -99.;
    v_bPartons_E.at(i)       = -99.;
    v_bPartons_M.at(i)       = -99.;
    v_bPartons_Origin.at(i)  = "EMPTY";

    v_SelectedJets_Lorentz.at(i) = TLorentzVector();
    v_SelectedJets_pt.at(i)      = -99.;
    v_SelectedJets_eta.at(i)     = -99.;
    v_SelectedJets_phi.at(i)     = -99.;
    v_SelectedJets_E.at(i)       = -99.;
    v_SelectedJets_M.at(i)       = -99.;
    v_SelectedJets_bTag.at(i)    = -99;

    v_bIndexMatchedToJet.at(i) = -99;
    v_JetMatchedOrigin.at(i)   = "EMPTY";
  }
  return;
}

//---------------------
// PRIVATE VIRTUAL Print information for each event
//---------------------
void RecoLevel::PrintEventInfo() {

  std::cout << "Event " << eventNum << std::endl;
  std::cout << "Passes cuts (true = 1): " << passSelection << std::endl;

  std::cout << " -----------" << std::endl;
  std::cout << " TRUTH LEVEL" << std::endl;
  std::cout << " -----------" << std::endl;

  if (m_sig) {
    std::cout << "  CHARGED HIGGS: (Charge,pt,eta,phi,E,M) = (" << truth_Hc_charge << ", " << truth_Hc_pt << ", ";
    std::cout << truth_Hc_eta << ", " << truth_Hc_phi << ", " << truth_Hc_E << ", " << truth_Hc_M << ")" << std::endl;
  }

  std::cout << "  TOP: (pt,eta,phi,E,M) = (" << truth_top_pt << ", " << truth_top_eta << ", ";
  std::cout << truth_top_phi << ", " << truth_top_E << ", " << truth_top_M << ")" << std::endl;

  std::cout << "  TBAR: (pt,eta,phi,E,M) = (" << truth_tbar_pt << ", " << truth_tbar_eta << ", ";
  std::cout << truth_tbar_phi << ", " << truth_tbar_E << ", " << truth_tbar_M << ")" << std::endl;

  std::cout << "  POS LEP: (pt,eta,phi,E,M) = (" << truth_poslep_pt << ", " << truth_poslep_eta << ", ";
  std::cout << truth_poslep_phi << ", " << truth_poslep_E << ", " << truth_poslep_M << ")" << std::endl;

  std::cout << "  NEG LEP: (pt,eta,phi,E,M) = (" << truth_neglep_pt << ", " << truth_neglep_eta << ", ";
  std::cout << truth_neglep_phi << ", " << truth_neglep_E << ", " << truth_neglep_M << ")" << std::endl;

  std::cout << "  NU: (pt,eta,phi,E,M) = (" << truth_nu_pt << ", " << truth_nu_eta << ", ";
  std::cout << truth_nu_phi << ", " << truth_nu_E << ", " << truth_nu_M << ")" << std::endl;

  std::cout << "  NUBAR: (pt,eta,phi,E,M) = (" << truth_nubar_pt << ", " << truth_nubar_eta << ", ";
  std::cout << truth_nubar_phi << ", " << truth_nubar_E << ", " << truth_nubar_M << ")" << std::endl;

  for (int i = 0; i < v_bPartons_pt.size(); i++) {
    std::cout << "  BQUARK " << i << ": (pt,eta,phi,E,M,origin) = (" << v_bPartons_pt[i] << ", " << v_bPartons_eta[i] << ", ";
    std::cout << v_bPartons_phi[i] << ", " << v_bPartons_E[i] << ", " << v_bPartons_M[i] << ", " << v_bPartons_Origin[i] << ")" << std::endl;
  }

  std::cout << " -----------" << std::endl;
  std::cout << " RECO LEVEL" << std::endl;
  std::cout << " -----------" << std::endl;

  std::cout << "  POS LEP: (pt,eta,phi,E,M) = (" << delphes_poslep_pt << ", " << delphes_poslep_eta << ", ";
  std::cout << delphes_poslep_phi << ", " << delphes_poslep_E << ", " << delphes_poslep_M << ")" << std::endl;

  std::cout << "  NEG LEP: (pt,eta,phi,E,M) = (" << delphes_neglep_pt << ", " << delphes_neglep_eta << ", ";
  std::cout << delphes_neglep_phi << ", " << delphes_neglep_E << ", " << delphes_neglep_M << ")" << std::endl;

  for (int j = 0; j < v_SelectedJets_pt.size(); j++) {
    std::cout << "  BJET " << j << ": (pt,eta,phi,E,M,btag,IndexOfMatchedQuark,OriginOfMatchQuark) = (";
    std::cout << v_SelectedJets_pt[j] << ", " << v_SelectedJets_eta[j] << ", " << v_SelectedJets_phi[j] << ", " << v_SelectedJets_E[j] << ", ";
    std::cout << v_SelectedJets_M[j] << v_SelectedJets_bTag[j] << ", " << v_bIndexMatchedToJet[j] << ", " << v_JetMatchedOrigin[j] << ")" << std::endl;
  }

  return;
}

//---------------------
// PUBLIC VIRTUAL Analyse event
//---------------------
void RecoLevel::Event(const char* inputFilename, const char* outputFilename, int evtNum, int numberOfBTags) {

  eventNum = evtNum;
  numJets = v_AllJetsLorentz.size();
  numBTags = numberOfBTags;

  InitialiseEventVariables();

  std::vector<TLorentzVector>::iterator allJets_itr = v_AllJetsLorentz.begin();
  std::vector<TLorentzVector>::iterator allJets_end = v_AllJetsLorentz.end();
  for ( ; allJets_itr != allJets_end; ++allJets_itr) {
    int i_jet = std::distance(v_AllJetsLorentz.begin(), allJets_itr);
    v_AllJets_pt.at(i_jet)   = (*allJets_itr).Pt();
    v_AllJets_eta.at(i_jet)  = (*allJets_itr).Eta();
    v_AllJets_phi.at(i_jet)  = (*allJets_itr).Phi();
    v_AllJets_E.at(i_jet)    = (*allJets_itr).E();
    v_AllJets_M.at(i_jet)    = (*allJets_itr).M();
  }

  std::vector<Jet*>::iterator allJ_itr = v_AllJets.begin();
  std::vector<Jet*>::iterator allJ_end = v_AllJets.end();
  for ( ; allJ_itr != allJ_end; ++allJ_itr) {
    int i_jet = std::distance(v_AllJets.begin(), allJ_itr);
    v_AllJets_bTag.at(i_jet) = (*allJ_itr)->BTag;
  }

  if (m_MissingET != NULL && v_SelectedJets.size() == 4 && v_RecoLeptonsLorentz.size() == 2) {
    passSelection = 1; // event passes cuts
    n_EventsPassingCuts++;

    //--------------
    // Fill truth level variables
    //--------------
    SetTruthVariables();

    //----------------
    // Set reco selected jets
    //----------------
    std::vector<std::pair<int, TLorentzVector> >::iterator jet_itr = v_SelectedJetsLorentz.begin();
    std::vector<std::pair<int, TLorentzVector> >::iterator jet_end = v_SelectedJetsLorentz.end();
    for ( ; jet_itr != jet_end; ++jet_itr ) {
      int i_jet = std::distance(v_SelectedJetsLorentz.begin(), jet_itr);
      v_SelectedJets_Lorentz.at(i_jet) = jet_itr->second;
      v_SelectedJets_pt.at(i_jet)      = v_SelectedJets_Lorentz.at(i_jet).Pt();
      v_SelectedJets_eta.at(i_jet)     = v_SelectedJets_Lorentz.at(i_jet).Eta();
      v_SelectedJets_phi.at(i_jet)     = v_SelectedJets_Lorentz.at(i_jet).Phi();
      v_SelectedJets_E.at(i_jet)       = v_SelectedJets_Lorentz.at(i_jet).E();
      v_SelectedJets_M.at(i_jet)       = v_SelectedJets_Lorentz.at(i_jet).M();
    }

    //---------------
    // Match selected reco jets to parton-level b-quarks
    //---------------
    MatchJetsToQuarks();
    numMatched = v_JetQuarkPair.size();

    for (int i_jm = 0; i_jm < v_JetQuarkPair.size(); i_jm++) {
      v_bIndexMatchedToJet.at((v_JetQuarkPair.at(i_jm)).first) = (v_JetQuarkPair.at(i_jm)).second;
    }
    for (int i_j = 0; i_j < v_SelectedJets.size(); i_j++) {
      v_SelectedJets_bTag.at(i_j) = (v_SelectedJets.at(i_j)).second->BTag;
      if (v_bIndexMatchedToJet.at(i_j) != -99 ) v_JetMatchedOrigin.at(i_j) = v_bPartons_Origin.at(v_bIndexMatchedToJet.at(i_j));
    }

    //----------------
    // Set reco charged leptons
    //----------------

    if (v_RecoLeptonsLorentz.size() == 2) {
      if (v_RecoLeptonsLorentz[0].first.second == +1) {
        delphes_poslep = v_RecoLeptonsLorentz[0].second;
        delphes_neglep = v_RecoLeptonsLorentz[1].second;
      }
      else if (v_RecoLeptonsLorentz[1].first.second == +1) {
        delphes_poslep = v_RecoLeptonsLorentz[1].second;
        delphes_neglep = v_RecoLeptonsLorentz[0].second;
      }
    }

    delphes_poslep_pt  = delphes_poslep.Pt();
    delphes_poslep_eta = delphes_poslep.Eta();
    delphes_poslep_phi = delphes_poslep.Phi();
    delphes_poslep_E   = delphes_poslep.E();
    delphes_poslep_M   = delphes_poslep.M();

    delphes_neglep_pt  = delphes_neglep.Pt();
    delphes_neglep_eta = delphes_neglep.Eta();
    delphes_neglep_phi = delphes_neglep.Phi();
    delphes_neglep_E   = delphes_neglep.E();
    delphes_neglep_M   = delphes_neglep.M();

    //----------------
    // Set ETmiss
    //----------------
    delphes_etmiss = m_MissingET->P4();
    delphes_etmiss_met = m_MissingET->MET;
    delphes_etmiss_eta = m_MissingET->Eta;
    delphes_etmiss_phi = m_MissingET->Phi;

  }

  int n_sameOrder = 0;
  if(passSelection == 1) {
    for (int i = 0; i < 4; i++) {
      if (v_bIndexMatchedToJet.at(i) != -1) {
        if (i == v_bIndexMatchedToJet.at(i)) n_sameOrder++;
      }
    }
  }

  if (numMatched == 3) n_EventsWith3Matches++;
  else if (numMatched == 4) n_EventsWith4Matches++;

  if (m_print) PrintEventInfo();

  //fill tree
  tree->Fill();

  v_JetQuarkPair.clear();

  return;
}

void RecoLevel::PrintInfo() {
  std::cout << "Number of events passing cuts: " << n_EventsPassingCuts*m_norm << std::endl;
  std::cout << "Number of events with 3 or 4 matches: respectively " << n_EventsWith3Matches*m_norm << ", " << n_EventsWith4Matches*m_norm << std::endl;
  return;
}
