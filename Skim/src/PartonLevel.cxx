#include "include/PartonLevel.h"

bool SortByDeltaR(std::pair<std::string,double> obj1, std::pair<std::string,double> obj2) {
  return (obj1.second > obj2.second);
}

//---------------------
// PUBLIC VIRTUAL Set tree branches reading for filling
//---------------------
void PartonLevel::BookTreeBranches() {
  std::cout << "Creating PartonResults tree..." << std::endl;

  tree = new TTree("PartonResults","PartonResults");
  tree->Branch("EventNumber", &eventNum, "EventNumber/I");

  if (m_sig) {
    tree->Branch("Truth_Hc_Lorentz", "TLorentzVector", &truth_Hc);
    tree->Branch("Truth_Hc_Charge", &truth_Hc_charge, "Truth_Hc_charge/I");
    tree->Branch("Truth_Hc_pt", &truth_Hc_pt, "Truth_Hc_pt/F");
    tree->Branch("Truth_Hc_eta", &truth_Hc_eta, "Truth_Hc_eta/F");
    tree->Branch("Truth_Hc_phi", &truth_Hc_phi, "Truth_Hc_phi/F");
    tree->Branch("Truth_Hc_E", &truth_Hc_E, "Truth_Hc_E/F");
    tree->Branch("Truth_Hc_M", &truth_Hc_M, "Truth_Hc_M/F");
  }

  tree->Branch("Truth_top_Lorentz", "TLorentzVector", &truth_top);
  tree->Branch("Truth_top_pt", &truth_top_pt, "Truth_top_pt/F");
  tree->Branch("Truth_top_eta", &truth_top_eta, "Truth_top_eta/F");
  tree->Branch("Truth_top_phi", &truth_top_phi, "Truth_top_phi/F");
  tree->Branch("Truth_top_E", &truth_top_E, "Truth_top_E/F");
  tree->Branch("Truth_top_M", &truth_top_M, "Truth_top_M/F");

  tree->Branch("Truth_tbar_Lorentz", "TLorentzVector", &truth_tbar);
  tree->Branch("Truth_tbar_pt", &truth_tbar_pt, "Truth_tbar_pt/F");
  tree->Branch("Truth_tbar_eta", &truth_tbar_eta, "Truth_tbar_eta/F");
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
  tree->Branch("Truth_nu_eta", &truth_nu_eta, "Truth_nu_eta/F");
  tree->Branch("Truth_nu_phi", &truth_nu_phi, "Truth_nu_phi/F");
  tree->Branch("Truth_nu_E", &truth_nu_E, "Truth_nu_E/F");
  tree->Branch("Truth_nu_M", &truth_nu_M, "Truth_nu_M/F");

  tree->Branch("Truth_nubar_Lorentz", "TLorentzVector", &truth_nubar);
  tree->Branch("Truth_nubar_pt", &truth_nubar_pt, "Truth_nubar_pt/F");
  tree->Branch("Truth_nubar_eta", &truth_nubar_eta, "Truth_nubar_eta/F");
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

  if (m_sig) {
    tree->Branch("Truth_bFromHiggs_Lorentz", "TLorentzVector", &bFromHiggs);
    tree->Branch("Truth_bFromHiggs_Index", &bFromHiggs_index, "Truth_bFromHiggs_Index/I");
    tree->Branch("Truth_bFromHiggs_pt", &bFromHiggs_pt, "Truth_bFromHiggs_pt/F");
    tree->Branch("Truth_bFromHiggs_eta",  &bFromHiggs_eta, "Truth_bFromHiggs_eta/F");
    tree->Branch("Truth_bFromHiggs_phi", &bFromHiggs_phi, "Truth_bFromHiggs_phi/F");
    tree->Branch("Truth_bFromHiggs_E", &bFromHiggs_E, "Truth_bFromHiggs_E/F");
    tree->Branch("Truth_bFromHiggs_M", &bFromHiggs_M, "Truth_bFromHiggs_M/F");

    tree->Branch("Truth_bFromTopHiggs_Lorentz", "TLorentzVector", &bFromTopHiggs);
    tree->Branch("Truth_bFromTopHiggs_Index", &bFromTopHiggs_index, "Truth_bFromTopHiggs_Index/I");
    tree->Branch("Truth_bFromTopHiggs_pt", &bFromTopHiggs_pt, "Truth_bFromTopHiggs_pt/F");
    tree->Branch("Truth_bFromTopHiggs_eta",  &bFromTopHiggs_eta, "Truth_bFromTopHiggs_eta/F");
    tree->Branch("Truth_bFromTopHiggs_phi", &bFromTopHiggs_phi, "Truth_bFromTopHiggs_phi/F");
    tree->Branch("Truth_bFromTopHiggs_E", &bFromTopHiggs_E, "Truth_bFromTopHiggs_E/F");
    tree->Branch("Truth_bFromTopHiggs_M", &bFromTopHiggs_M, "Truth_bFromTopHiggs_M/F");

    tree->Branch("Truth_bFromTop_Lorentz", "TLorentzVector", &bFromTop);
    tree->Branch("Truth_bFromTop_Index", &bFromTop_index, "Truth_bFromTop_Index/I");
    tree->Branch("Truth_bFromTop_pt", &bFromTop_pt, "Truth_bFromTop_pt/F");
    tree->Branch("Truth_bFromTop_eta",  &bFromTop_eta, "Truth_bFromTop_eta/F");
    tree->Branch("Truth_bFromTop_phi", &bFromTop_phi, "Truth_bFromTop_phi/F");
    tree->Branch("Truth_bFromTop_E", &bFromTop_E, "Truth_bFromTop_E/F");
    tree->Branch("Truth_bFromTop_M", &bFromTop_M, "Truth_bFromTop_M/F");

    tree->Branch("Truth_bFromGluon_Lorentz", "TLorentzVector", &bFromGluon);
    tree->Branch("Truth_bFromGluon_Index", &bFromGluon_index, "Truth_bFromGluon_Index/I");
    tree->Branch("Truth_bFromGluon_pt", &bFromGluon_pt, "Truth_bFromGluon_pt/F");
    tree->Branch("Truth_bFromGluon_eta",  &bFromGluon_eta, "Truth_bFromGluon_eta/F");
    tree->Branch("Truth_bFromGluon_phi", &bFromGluon_phi, "Truth_bFromGluon_phi/F");
    tree->Branch("Truth_bFromGluon_E", &bFromGluon_E, "Truth_bFromGluon_E/F");
    tree->Branch("Truth_bFromGluon_M", &bFromGluon_M, "Truth_bFromGluon_M/F");
  }
  else {
    tree->Branch("Truth_bFromTop1_Lorentz", "TLorentzVector", &bFromTop1);
    tree->Branch("Truth_bFromTop1_Index", &bFromTop1_index, "Truth_bFromTop1_Index/I");
    tree->Branch("Truth_bFromTop1_pt", &bFromTop1_pt, "Truth_bFromTop1_pt/F");
    tree->Branch("Truth_bFromTop1_eta",  &bFromTop1_eta, "Truth_bFromTop1_eta/F");
    tree->Branch("Truth_bFromTop1_phi", &bFromTop1_phi, "Truth_bFromTop1_phi/F");
    tree->Branch("Truth_bFromTop1_E", &bFromTop1_E, "Truth_bFromTop1_E/F");
    tree->Branch("Truth_bFromTop1_M", &bFromTop1_M, "Truth_bFromTop1_M/F");

    tree->Branch("Truth_bFromTop2_Lorentz", "TLorentzVector", &bFromTop2);
    tree->Branch("Truth_bFromTop2_Index", &bFromTop2_index, "Truth_bFromTop2_Index/I");
    tree->Branch("Truth_bFromTop2_pt", &bFromTop2_pt, "Truth_bFromTop2_pt/F");
    tree->Branch("Truth_bFromTop2_eta", &bFromTop2_eta, "Truth_bFromTop2_eta/F");
    tree->Branch("Truth_bFromTop2_phi", &bFromTop2_phi, "Truth_bFromTop2_phi/F");
    tree->Branch("Truth_bFromTop2_E", &bFromTop2_E, "Truth_bFromTop2_E/F");
    tree->Branch("Truth_bFromTop2_M", &bFromTop2_M, "Truth_bFromTop2_M/F");

    tree->Branch("Truth_bFromGluon1_Lorentz", "TLorentzVector", &bFromGluon1);
    tree->Branch("Truth_bFromGluon1_Index", &bFromGluon1_index, "Truth_bFromGluon1_Index/I");
    tree->Branch("Truth_bFromGluon1_pt", &bFromGluon1_pt, "Truth_bFromGluon1_pt/F");
    tree->Branch("Truth_bFromGluon1_eta",  &bFromGluon1_eta, "Truth_bFromGluon1_eta/F");
    tree->Branch("Truth_bFromGluon1_phi", &bFromGluon1_phi, "Truth_bFromGluon1_phi/F");
    tree->Branch("Truth_bFromGluon1_E", &bFromGluon1_E, "Truth_bFromGluon1_E/F");
    tree->Branch("Truth_bFromGluon1_M", &bFromGluon1_M, "Truth_bFromGluon1_M/F");

    tree->Branch("Truth_bFromGluon2_Lorentz", "TLorentzVector", &bFromGluon2);
    tree->Branch("Truth_bFromGluon2_Index", &bFromGluon2_index, "Truth_bFromGluon2_Index/I");
    tree->Branch("Truth_bFromGluon2_pt", &bFromGluon2_pt, "Truth_bFromGluon2_pt/F");
    tree->Branch("Truth_bFromGluon2_eta",  &bFromGluon2_eta, "Truth_bFromGluon2_eta/F");
    tree->Branch("Truth_bFromGluon2_phi", &bFromGluon2_phi, "Truth_bFromGluon2_phi/F");
    tree->Branch("Truth_bFromGluon2_E", &bFromGluon2_E, "Truth_bFromGluon2_E/F");
    tree->Branch("Truth_bFromGluon2_M", &bFromGluon2_M, "Truth_bFromGluon2_M/F");
  }

  //---------------
  // DeltaR between pairs of b-partons
  //---------------
  if (m_sig) {
    tree->Branch("DeltaR_bH_bt", &deltaR_bH_bt, "DeltaR_bH_bt/F");
    tree->Branch("DeltaR_bH_bHt", &deltaR_bH_bHt, "DeltaR_bH_bHt/F");
    tree->Branch("DeltaR_bH_bg", &deltaR_bH_bg, "DeltaR_bH_bg/F");
    tree->Branch("DeltaR_bHt_bt", &deltaR_bHt_bt, "DeltaR_bHt_bt/F");
    tree->Branch("DeltaR_bHt_bg", &deltaR_bHt_bg, "DeltaR_bHt_bg/F");
    tree->Branch("DeltaR_bt_bg", &deltaR_bt_bg, "DeltaR_bt_bg/F");
  }
  else {
    tree->Branch("deltaR_bt1_bt2", &deltaR_bt1_bt2, "DeltaR_bt1_bt2/F");
    tree->Branch("DeltaR_bt1_bg1", &deltaR_bt1_bg1, "DeltaR_bt1_bg1/F");
    tree->Branch("DeltaR_bt2_bg2", &deltaR_bt1_bg2, "DeltaR_bt1_bg2/F");
    tree->Branch("DeltaR_bt2_bg1", &deltaR_bt2_bg1, "DeltaR_bt2_bg1/F");
    tree->Branch("DeltaR_bt2_bg2", &deltaR_bt2_bg2, "DeltaR_bt2_bg2/F");
    tree->Branch("DeltaR_bg1_bg2", &deltaR_bg1_bg2, "DeltaR_bg1_bg2/F");
  }
  // gives the origin of the two b-partons in each case
  tree->Branch("Identifier_ClosestPair_bb", &closestPair);
  tree->Branch("Identifier_2ndClosestPair_bb", &secondPair);
  tree->Branch("Identifier_3rdClosestPair_bb", &thirdPair);
  tree->Branch("Identifier_4thClosestPair_bb", &fourthPair);
  tree->Branch("Identifier_5thClosestPair_bb", &fifthPair);
  tree->Branch("Identifier_6thClosestPair_bb", &sixthPair);

  return;
}

//---------------------
// PRIVATE VIRTUAL Set tree variables to default value
//---------------------
void PartonLevel::InitialiseEventVariables() {

  for (int i = 0; i < 4; i++) v_bPartons_Lorentz.at(i) = TLorentzVector();

  deltaR_bH_bHt = -99.;
  deltaR_bH_bt  = -99.;
  deltaR_bH_bg  = -99.;
  deltaR_bHt_bt = -99.;
  deltaR_bHt_bg = -99.;
  deltaR_bt_bg  = -99.;

  deltaR_bt1_bt2 = -99.;
  deltaR_bt1_bg1  = -99.;
  deltaR_bt1_bg2  = -99.;
  deltaR_bt2_bg1 = -99.;
  deltaR_bt2_bg2 = -99.;
  deltaR_bg1_bg2  = -99.;

  closestPair = "EMPTY";
  secondPair  = "EMPTY";
  thirdPair   = "EMPTY";
  fourthPair  = "EMPTY";
  fifthPair   = "EMPTY";
  sixthPair   = "EMPTY";

  truth_Hc        = TLorentzVector();
  truth_top       = TLorentzVector();
  truth_tbar      = TLorentzVector();
  truth_nu        = TLorentzVector();
  truth_nubar     = TLorentzVector();
  truth_poslep    = TLorentzVector();
  truth_neglep    = TLorentzVector();
  bFromTop        = TLorentzVector();
  bFromTopHiggs   = TLorentzVector();
  bFromHiggs      = TLorentzVector();
  bFromGluon      = TLorentzVector();
  bFromTop_index      = -1;
  bFromTopHiggs_index = -1;
  bFromHiggs_index    = -1;
  bFromGluon_index    = -1;

  bFromTop1     = TLorentzVector();
  bFromTop2     = TLorentzVector();
  bFromGluon1   = TLorentzVector();
  bFromGluon2   = TLorentzVector();
  bFromTop1_index   = -1;
  bFromTop2_index   = -1;
  bFromGluon1_index = -1;
  bFromGluon2_index = -1;

  return;
}

//---------------------
// PRIVATE VIRTUAL Print information for each event
//---------------------
void PartonLevel::PrintEventInfo() {

  std::cout << "Event " << eventNum << std::endl;

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
    std::cout << "  BQUARK " << i << ": (pt,eta,phi,E,M) = (" << v_bPartons_pt[i] << ", " << v_bPartons_eta[i] << ", ";
    std::cout << v_bPartons_phi[i] << ", " << v_bPartons_E[i] << ", " << v_bPartons_M[i] << ")" << std::endl;
  }

  if (m_sig) {
    std::cout << "  b FROM HIGGS: (pt,eta,phi,E,M,index) = (" << bFromHiggs_pt << ", " << bFromHiggs_eta << ", ";
    std::cout << bFromHiggs_phi << ", " << bFromHiggs_E << ", " << bFromHiggs_M << ", " << bFromHiggs_index << ")" << std::endl;

    std::cout << "  b FROM TOP AND HIGGS: (pt,eta,phi,E,M) = (" << bFromTopHiggs_pt << ", " << bFromTopHiggs_eta << ", ";
    std::cout << bFromTopHiggs_phi << ", " << bFromTopHiggs_E << ", " << bFromTopHiggs_M << ", " << bFromTopHiggs_index << ")" << std::endl;

    std::cout << "  b FROM TOP: (pt,eta,phi,E,M) = (" << bFromTop_pt << ", " << bFromTop_eta << ", ";
    std::cout << bFromTop_phi << ", " << bFromTop_E << ", " << bFromTop_M << ", " << bFromTop_index << ")" << std::endl;

    std::cout << "  b FROM GLUON: (pt,eta,phi,E,M) = (" << bFromGluon_pt << ", " << bFromGluon_eta << ", ";
    std::cout << bFromGluon_phi << ", " << bFromGluon_E << ", " << bFromGluon_M << ", " << bFromGluon_index << ")" << std::endl;
  }
  else {
    std::cout << "  b FROM TOP 1: (pt,eta,phi,E,M) = (" << bFromTop1_pt << ", " << bFromTop1_eta << ", ";
    std::cout << bFromTop1_phi << ", " << bFromTop1_E << ", " << bFromTop1_M << ", " << bFromTop1_index << ")" << std::endl;

    std::cout << "  b FROM TOP 2: (pt,eta,phi,E,M) = (" << bFromTop2_pt << ", " << bFromTop2_eta << ", ";
    std::cout << bFromTop2_phi << ", " << bFromTop2_E << ", " << bFromTop2_M << ", " << bFromTop2_index << ")" << std::endl;

    std::cout << "  b FROM GLUON 1: (pt,eta,phi,E,M) = (" << bFromGluon1_pt << ", " << bFromGluon1_eta << ", ";
    std::cout << bFromGluon1_phi << ", " << bFromGluon1_E << ", " << bFromGluon1_M << ", " << bFromGluon1_index << ")" << std::endl;

    std::cout << "  b FROM GLUON 2: (pt,eta,phi,E,M) = (" << bFromGluon2_pt << ", " << bFromGluon2_eta << ", ";
    std::cout << bFromGluon2_phi << ", " << bFromGluon2_E << ", " << bFromGluon2_M << ", " << bFromGluon2_index << ")" << std::endl;
  }

  std::cout << " -----------" << std::endl;
  std::cout << " TRUTH LEVEL DELTA R" << std::endl;
  std::cout << " -----------" << std::endl;
  if (m_sig) {
    std::cout << "  DR(bH,bt) = " << deltaR_bH_bt << " DR(bH,bHt) = " << deltaR_bH_bHt << " DR(bH,bg) = " << deltaR_bH_bg;
    std::cout << "  DR(bHt,bt) = " << deltaR_bHt_bt << " DR(bHt,bg) = " << deltaR_bHt_bg << " DR(bt,bg) = " << deltaR_bt_bg << std::endl;
  }
  else {
    std::cout << "  DR(bt1,bt2) = " << deltaR_bt1_bt2 << " DR(bt1,bg1) = " << deltaR_bt1_bg1 << " DR(bt1,bg2) = " << deltaR_bt1_bg2;
    std::cout << "  DR(bt2,bg1) = " << deltaR_bt2_bg1 << " DR(bt2,bg2) = " << deltaR_bt2_bg2 << " DR(bg1,bg2) = " << deltaR_bg1_bg2 << std::endl;
  }
  std::cout << "  Closest = " << closestPair << " 2nd = " << secondPair << " 3rd = " << thirdPair;
  std::cout << "  4th = " << fourthPair << " 5th = " << fifthPair << " 6th = " << sixthPair << std::endl;

  return;
}

//---------------------
// PUBLIC VIRTUAL Analyse event
//---------------------
void PartonLevel::Event(const char* inputFilename, const char* outputFilename, int evtNum, int numberOfBTags) {

  eventNum = evtNum;
  InitialiseEventVariables();

  //------------
  // TRUTH VARIABLES
  //------------
  SetTruthVariables();

  // signal sample
  if (m_sig) {

    bFromHiggs    = v_GenParticlesLorentz["bFromHiggs"];
    bFromTopHiggs = v_GenParticlesLorentz["bFromTopAndHiggs"];
    bFromTop      = v_GenParticlesLorentz["bFromTop"];
    bFromGluon    = v_GenParticlesLorentz["bFromGluon"];

    bFromTopHiggs_pt  = bFromTopHiggs.Pt();
    bFromTopHiggs_eta = bFromTopHiggs.Eta();
    bFromTopHiggs_phi = bFromTopHiggs.Phi();
    bFromTopHiggs_E   = bFromTopHiggs.E();
    bFromTopHiggs_M   = bFromTopHiggs.M();

    bFromTop_pt  = bFromTop.Pt();
    bFromTop_eta = bFromTop.Eta();
    bFromTop_phi = bFromTop.Phi();
    bFromTop_E   = bFromTop.E();
    bFromTop_M   = bFromTop.M();

    bFromHiggs_pt  = bFromHiggs.Pt();
    bFromHiggs_eta = bFromHiggs.Eta();
    bFromHiggs_phi = bFromHiggs.Phi();
    bFromHiggs_E   = bFromHiggs.E();
    bFromHiggs_M   = bFromHiggs.M();

    bFromGluon_pt  = bFromGluon.Pt();
    bFromGluon_eta = bFromGluon.Eta();
    bFromGluon_phi = bFromGluon.Phi();
    bFromGluon_E   = bFromGluon.E();
    bFromGluon_M   = bFromGluon.M();

    std::vector<std::pair<std::string,TLorentzVector> >::iterator b_itr = v_bQuarksLorentz.begin();
    std::vector<std::pair<std::string,TLorentzVector> >::iterator b_end = v_bQuarksLorentz.end();
    for ( ; b_itr != b_end; ++b_itr) {
      int i_b = std::distance(v_bQuarksLorentz.begin(), b_itr);
      if ((b_itr->first).compare("bFromHiggs") == 0) bFromHiggs_index = i_b;
      if ((b_itr->first).compare("bFromTopAndHiggs") == 0) bFromTopHiggs_index = i_b;
      if ((b_itr->first).compare("bFromTop") == 0) bFromTop_index = i_b;
      if ((b_itr->first).compare("bFromGluon") == 0) bFromGluon_index = i_b;
    }

    //------------------
    // Calculate DeltaR between all pairs of b-quarks and determine what is closest
    //------------------
    deltaR_bH_bHt = bFromHiggs.DeltaR(bFromTopHiggs);
    deltaR_bH_bt  = bFromHiggs.DeltaR(bFromTop);
    deltaR_bH_bg  = bFromHiggs.DeltaR(bFromGluon);
    deltaR_bHt_bt = bFromTopHiggs.DeltaR(bFromTop);
    deltaR_bHt_bg = bFromTopHiggs.DeltaR(bFromGluon);
    deltaR_bt_bg  = bFromTop.DeltaR(bFromGluon);

    deltaR_bquarks.at(0) = std::make_pair("#DeltaR(b_{H},b_{Ht})", deltaR_bH_bHt);
    deltaR_bquarks.at(1) = std::make_pair("#DeltaR(b_{H},b_{t})", deltaR_bH_bt);
    deltaR_bquarks.at(2) = std::make_pair("#DeltaR(b_{H},b_{g})", deltaR_bH_bg);
    deltaR_bquarks.at(3) = std::make_pair("#DeltaR(b_{Ht},b_{t})", deltaR_bHt_bt);
    deltaR_bquarks.at(4) = std::make_pair("#DeltaR(b_{Ht},b_{g})", deltaR_bHt_bg);
    deltaR_bquarks.at(5) = std::make_pair("#DeltaR(b_{t},b_{g})", deltaR_bt_bg);
  }
  // background sample
  else {
    bFromTop1   = v_GenParticlesLorentz["bFromTop1"];
    bFromTop2   = v_GenParticlesLorentz["bFromTop2"];
    bFromGluon1 = v_GenParticlesLorentz["bFromGluon1"];
    bFromGluon2 = v_GenParticlesLorentz["bFromGluon2"];

    bFromTop1_pt  = bFromTop1.Pt();
    bFromTop1_eta = bFromTop1.Eta();
    bFromTop1_phi = bFromTop1.Phi();
    bFromTop1_E   = bFromTop1.E();
    bFromTop1_M   = bFromTop1.M();

    bFromTop2_pt  = bFromTop2.Pt();
    bFromTop2_eta = bFromTop2.Eta();
    bFromTop2_phi = bFromTop2.Phi();
    bFromTop2_E   = bFromTop2.E();
    bFromTop2_M   = bFromTop2.M();

    bFromGluon1_pt  = bFromGluon1.Pt();
    bFromGluon1_eta = bFromGluon1.Eta();
    bFromGluon1_phi = bFromGluon1.Phi();
    bFromGluon1_E   = bFromGluon1.E();
    bFromGluon1_M   = bFromGluon1.M();

    bFromGluon2_pt  = bFromGluon2.Pt();
    bFromGluon2_eta = bFromGluon2.Eta();
    bFromGluon2_phi = bFromGluon2.Phi();
    bFromGluon2_E   = bFromGluon2.E();
    bFromGluon2_M   = bFromGluon2.M();

    std::vector<std::pair<std::string,TLorentzVector> >::iterator b_itr = v_bQuarksLorentz.begin();
    std::vector<std::pair<std::string,TLorentzVector> >::iterator b_end = v_bQuarksLorentz.end();
    for ( ; b_itr != b_end; ++b_itr) {
      int i_b = std::distance(v_bQuarksLorentz.begin(), b_itr);
      if ((b_itr->first).compare("bFromTop1") == 0) bFromTop1_index = i_b;
      if ((b_itr->first).compare("bFromTop2") == 0) bFromTop2_index = i_b;
      if ((b_itr->first).compare("bFromGluon1") == 0) bFromGluon1_index = i_b;
      if ((b_itr->first).compare("bFromGluon2") == 0) bFromGluon2_index = i_b;
    }

    //------------------
    // Calculate DeltaR between all pairs of b-quarks and determine what is closest
    //------------------
    deltaR_bt1_bt2 = bFromTop1.DeltaR(bFromTop2);
    deltaR_bt1_bg1 = bFromTop1.DeltaR(bFromGluon1);
    deltaR_bt1_bg2 = bFromTop1.DeltaR(bFromGluon2);
    deltaR_bt2_bg1 = bFromTop2.DeltaR(bFromGluon1);
    deltaR_bt2_bg2 = bFromTop2.DeltaR(bFromGluon2);
    deltaR_bg1_bg2 = bFromGluon1.DeltaR(bFromGluon2);

    deltaR_bquarks.at(0) = std::make_pair("#DeltaR(b_{t_{1}},b_{t_{2}})", deltaR_bt1_bt2);
    deltaR_bquarks.at(1) = std::make_pair("#DeltaR(b_{t_{1}},b_{g_{1}})", deltaR_bt1_bg1);
    deltaR_bquarks.at(2) = std::make_pair("#DeltaR(b_{t_{1}},b_{g_{2}})", deltaR_bt1_bg2);
    deltaR_bquarks.at(3) = std::make_pair("#DeltaR(b_{t_{2}},b_{g_{1}})", deltaR_bt2_bg1);
    deltaR_bquarks.at(4) = std::make_pair("#DeltaR(b_{t_{2}},b_{g_{2}})", deltaR_bt2_bg2);
    deltaR_bquarks.at(5) = std::make_pair("#DeltaR(b_{g_{1}},b_{g_{2}})", deltaR_bg1_bg2);
  }

  std::sort(deltaR_bquarks.begin(), deltaR_bquarks.end(), SortByDeltaR);

  closestPair = deltaR_bquarks[0].first;
  secondPair  = deltaR_bquarks[1].first;
  thirdPair   = deltaR_bquarks[2].first;
  fourthPair  = deltaR_bquarks[3].first;
  fifthPair   = deltaR_bquarks[4].first;
  sixthPair   = deltaR_bquarks[5].first;

  // check whether variables are being filled correctly in tree
  if (m_print) PrintEventInfo();

  tree->Fill();

  return;
}

void PartonLevel::PrintInfo() {
  std::cout << "All events processed." << std::endl;
}
