#ifndef LIMITSHISTMAKER_H
#define LIMITSHISTMAKER_H

#include "include/RootIncludes.h"

class LimitsHistMaker {

private:

  void BookHists();
  void GetBranches();
  void SetTMVAReader();
  void GetBDTResult();

  float GetMinv_b_ClosestTop(TLorentzVector* b);
  float GetMinv_JetsWithMinDR();
  float GetCosTheta();
  float GetHt();
  float GetMinvMin_lep_b(std::string charge);
  float GetDE_NonTopJets();
  float GetCentrality();
  float GetDR(TLorentzVector* a, TLorentzVector* b);
  float GetMinvPair(TLorentzVector* a, TLorentzVector* b);

  void WriteHists();
  void DeleteChainVariables();
  void DeleteHists();

  std::string m_inputFileName;
  std::string m_outputFileName;
  std::string m_typeBDT;
  std::string m_WeightsFile;

	TTree* m_inputTree;
	TTree* m_outputTree;

	TFile* m_inputFile;
	TFile* m_outputFile;

  TMVA::Reader* m_reader;

  bool m_haveRecoBDTWeightCut;
  float m_recoBDTWeightCutValue;

  bool m_haveBDTWeightCut;
  float m_BDTWeightCutValue;

  float m_norm;

  //-------------------
  // Histograms
  //-------------------
  TH1F* h_numJets;
  TH1F* h_maxRecoBDTWeight;
  // TH1F* h_minv_bLeading_ClosestTop;
  TH1F* h_minv_bH_ClosestTop;
  TH1F* h_minv_minDRJets;
  TH1F* h_Ht;
  TH1F* h_minMinv_poslep_b;
  TH1F* h_minMinv_neglep_b;
  // TH1F* h_pt_leadingJet;
  TH1F* h_pt_bH;
  TH1F* h_minv_Higgs;
  TH1F* h_dE_NonTopJets;
  TH1F* h_dR_Tops;
  TH1F* h_dR_bH_tH;
  TH1F* h_dR_bLeading_TopLeading;
  TH1F* h_centrality;
  TH1F* h_dR_tH_lt;
  TH1F* h_dR_lep_maxDR_t_bH;
  TH1F* h_costheta_lepton_jet;
  TH1F* h_eta_bH;
  TH1F* h_eta_btH;
  TH1F* h_dR_bH_btH;
  TH1F* h_dR_bg_bt;
  TH1F* h_dEta_bg_bt;
  TH1F* h_dPhi_btH_bH;
  TH1F* h_dPhi_bt_bH;
  TH1F* h_dPhi_bt_btH;
  TH1F* h_DR_bg_bH;
  TH1F* h_dphi_poslep_neglep;
  TH1F* h_dphi_ltH_bH;
  TH1F* h_DR_ltH_bH;
  TH1F* h_dphi_ltH_bt;
  TH1F* h_DR_ltH_bt;

  //-------------------
  // Event variables
  //-------------------
  int evt, passSelection, numJets, numBTags, numMatched;
  int hasRequiredSIGMatches, hasRequiredBKGMatches;
  float maxRecoBDTWeight;

  std::vector<int>* v_AllJets_bTag;
  std::vector<float>* v_AllJets_pt;
  std::vector<float>* v_AllJets_eta;
  std::vector<float>* v_AllJets_phi;
  std::vector<float>* v_AllJets_E;
  std::vector<float>* v_AllJets_M;

  std::vector<TLorentzVector*> v_SelectedJets_Lorentz;
  std::vector<int>* v_SelectedJets_bTag;
  std::vector<float>* v_SelectedJets_pt;
  std::vector<float>* v_SelectedJets_eta;
  std::vector<float>* v_SelectedJets_phi;
  std::vector<float>* v_SelectedJets_E;
  std::vector<float>* v_SelectedJets_M;

  std::vector<int>* v_bIndexMatchedToJet;
  std::vector<std::string>* v_JetMatchedOrigin;

  TLorentzVector* delphes_poslep;
  float delphes_poslep_pt, delphes_poslep_eta, delphes_poslep_phi, delphes_poslep_E, delphes_poslep_M;

  TLorentzVector* delphes_neglep;
  float delphes_neglep_pt, delphes_neglep_eta, delphes_neglep_phi, delphes_neglep_E, delphes_neglep_M;

  TLorentzVector* delphes_etmiss;
  float delphes_etmiss_met, delphes_etmiss_eta, delphes_etmiss_phi;

  float maxBDTWeight;
  int isReconstructed;
  float maxNeutrinoWeight;
  // tops
  TLorentzVector* maxNW_top;
  TLorentzVector* maxNW_tbar;
  float maxNW_top_pt, maxNW_top_eta, maxNW_top_phi, maxNW_top_E, maxNW_top_M;
  float maxNW_tbar_pt, maxNW_tbar_eta, maxNW_tbar_phi, maxNW_tbar_E, maxNW_tbar_M;
  // bjets
  TLorentzVector* maxNW_b;
  TLorentzVector* maxNW_bbar;
  int maxNW_b_index, maxNW_bbar_index, maxNW_bH_index, maxNW_bg_index, maxNW_btH_index, maxNW_bt_index;
  float maxNW_b_pt, maxNW_b_eta, maxNW_b_phi, maxNW_b_E, maxNW_b_M;
  float maxNW_bbar_pt, maxNW_bbar_eta, maxNW_bbar_phi, maxNW_bbar_E, maxNW_bbar_M;
  // neutrinos
  TLorentzVector* maxNW_nu;
  TLorentzVector* maxNW_nubar;
  float maxNW_nu_pt, maxNW_nu_eta, maxNW_nu_phi, maxNW_nu_E, maxNW_nu_M;
  float maxNW_nubar_pt, maxNW_nubar_eta, maxNW_nubar_phi, maxNW_nubar_E, maxNW_nubar_M;

  TLorentzVector* topFromHiggs;
  TLorentzVector* nonHiggsLepton;
  TLorentzVector* higgsLepton;
  TLorentzVector* bFromHiggs;
  TLorentzVector* bFromTopAndHiggs;
  TLorentzVector* bFromTop;
  TLorentzVector* bFromGluon;

  //-------------------
  // Training sig/bkg variables
  //-------------------
  float m_maxRecoBDTWeight;
  float m_njets;
  int m_reconstructedCharge;
  float m_minv_bleading_closesttop;
  float m_minv_bH_closesttop;
  float m_minv_minDRJets;
  float m_ht;
  float m_minMinv_poslep_b;
  float m_minMinv_neglep_b;
  int m_nBTags;
  float m_leadingJetPt;
  float m_bHPt;
  float m_minv_Higgs;
  float m_DE_nonTopJets;
  float m_DR_tops;
  float m_DR_bH_tH;
  float m_DR_bleading_TopLeading;
  float m_centrality;
  float m_DR_tH_lt;
  float m_DR_lep_maxDR_t_bH;
  float m_costheta_lepton_jet;
  float m_eta_bH;
  float m_eta_btH;
  float m_deta_bg_bt;
  float m_dphi_btH_bH;
  float m_dphi_bt_bH;
  float m_DR_bH_btH;
  float m_DR_bg_bt;
  float m_dphi_bt_btH;
  float m_DR_bg_bH;
  float m_dphi_poslep_neglep;
  float m_dphi_ltH_bH;
  float m_DR_ltH_bH;
  float m_dphi_ltH_bt;
  float m_DR_ltH_bt;

  float m_BDTWeight;

public:

  LimitsHistMaker(std::string inputFileName, std::string outputFileName, std::string typeBDT, bool haveBDTWeightCut, float bdtWeightCutValue, float normalisation);
  ~LimitsHistMaker() {
    delete m_outputFile;
    delete m_inputFile;
  };

  void Run();

};

#endif
