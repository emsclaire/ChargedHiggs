#ifndef RECOTRAININGNTUPLEMAKER_H
#define RECOTRAININGNTUPLEMAKER_H

#include "include/RootIncludes.h"

class RecoTrainingNtupleMaker {

private:

  void ReadInputFiles();
  void GetBranches();
  void GetFunctionMap();
  void GetJetJetTrainingVariables(std::string function);
  void GetJetLepTrainingVariables(std::string function);
  void GetJetJetFunctions(std::string function);
  void GetJetLepFunctions(std::string function);
  float DPtPair(TLorentzVector* a, TLorentzVector* b);
  float DEtaPair(TLorentzVector* a, TLorentzVector* b);
  float DPhiPair(TLorentzVector* a, TLorentzVector* b);
  float DMPair(TLorentzVector* a, TLorentzVector* b);
  float DRPair(TLorentzVector* a, TLorentzVector* b);
  float PtPair(TLorentzVector* a, TLorentzVector* b);
  float MinvPair(TLorentzVector* a, TLorentzVector* b);
  float EtaPair(TLorentzVector* a, TLorentzVector* b);
  float PhiPair(TLorentzVector* a, TLorentzVector* b);
  float EPair(TLorentzVector* a, TLorentzVector* b);
  void DeleteChainVariables();

  std::string m_inputListFileName;
  std::string m_outputFileName;
  std::vector<std::string> v_inputFileList;

  TChain* delphes_chain;

  //-------------------
  // Delphes chain variables
  //-------------------
  int evt, passSelection, numJets, numBTags, numMatched;
  int hasRequiredSIGMatches, hasRequiredBKGMatches;

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

  //-------------------
  // Training ntuple variables
  //-------------------
  int m_evtNum;
  int m_nJets, m_nBTags;
  int m_correctmatch;

  int m_HcJet_index, m_tHcJet_index, m_tJet_index, m_gJet_index;
  float m_HcJet_pt, m_tHcJet_pt, m_tJet_pt, m_gJet_pt;
  float m_HcJet_eta, m_tHcJet_eta, m_tJet_eta, m_gJet_eta;
  float m_HcJet_phi, m_tHcJet_phi, m_tJet_phi, m_gJet_phi;
  float m_HcJet_e, m_tHcJet_e, m_tJet_e, m_gJet_e;
  float m_HcJet_tHcJet_dpt, m_HcJet_tJet_dpt, m_HcJet_gJet_dpt, m_tHcJet_tJet_dpt, m_tHcJet_gJet_dpt, m_tJet_gJet_dpt;
  float m_HcJet_tHcJet_deta, m_HcJet_tJet_deta, m_HcJet_gJet_deta, m_tHcJet_tJet_deta, m_tHcJet_gJet_deta, m_tJet_gJet_deta;
  float m_HcJet_tHcJet_dphi, m_HcJet_tJet_dphi, m_HcJet_gJet_dphi, m_tHcJet_tJet_dphi, m_tHcJet_gJet_dphi, m_tJet_gJet_dphi;
  float m_HcJet_tHcJet_dm, m_HcJet_tJet_dm, m_HcJet_gJet_dm, m_tHcJet_tJet_dm, m_tHcJet_gJet_dm, m_tJet_gJet_dm;
  float m_HcJet_tHcJet_dR, m_HcJet_tJet_dR, m_HcJet_gJet_dR, m_tHcJet_tJet_dR, m_tHcJet_gJet_dR, m_tJet_gJet_dR;
  float m_HcJet_tHcJet_ptpair, m_HcJet_tJet_ptpair, m_HcJet_gJet_ptpair, m_tHcJet_tJet_ptpair, m_tHcJet_gJet_ptpair, m_tJet_gJet_ptpair;
  float m_HcJet_tHcJet_etapair, m_HcJet_tJet_etapair, m_HcJet_gJet_etapair, m_tHcJet_tJet_etapair, m_tHcJet_gJet_etapair, m_tJet_gJet_etapair;
  float m_HcJet_tHcJet_phipair, m_HcJet_tJet_phipair, m_HcJet_gJet_phipair, m_tHcJet_tJet_phipair, m_tHcJet_gJet_phipair, m_tJet_gJet_phipair;
  float m_HcJet_tHcJet_epair, m_HcJet_tJet_epair, m_HcJet_gJet_epair, m_tHcJet_tJet_epair, m_tHcJet_gJet_epair, m_tJet_gJet_epair;
  float m_HcJet_PosLep_minv, m_tHcJet_PosLep_minv, m_tJet_PosLep_minv, m_gJet_PosLep_minv;
  float m_HcJet_NegLep_minv, m_tHcJet_NegLep_minv, m_tJet_NegLep_minv, m_gJet_NegLep_minv;
  float m_HcJet_PosLep_ptpair, m_tHcJet_PosLep_ptpair, m_tJet_PosLep_ptpair, m_gJet_PosLep_ptpair;
  float m_HcJet_NegLep_ptpair, m_tHcJet_NegLep_ptpair, m_tJet_NegLep_ptpair, m_gJet_NegLep_ptpair;
  float m_HcJet_PosLep_etapair, m_tHcJet_PosLep_etapair, m_tJet_PosLep_etapair, m_gJet_PosLep_etapair;
  float m_HcJet_NegLep_etapair, m_tHcJet_NegLep_etapair, m_tJet_NegLep_etapair, m_gJet_NegLep_etapair;
  float m_HcJet_PosLep_phipair, m_tHcJet_PosLep_phipair, m_tJet_PosLep_phipair, m_gJet_PosLep_phipair;
  float m_HcJet_NegLep_phipair, m_tHcJet_NegLep_phipair, m_tJet_NegLep_phipair, m_gJet_NegLep_phipair;
  float m_HcJet_PosLep_epair, m_tHcJet_PosLep_epair, m_tJet_PosLep_epair, m_gJet_PosLep_epair;
  float m_HcJet_NegLep_epair, m_tHcJet_NegLep_epair, m_tJet_NegLep_epair, m_gJet_NegLep_epair;

public:

  RecoTrainingNtupleMaker(std::string inputListFileName, std::string outputFileName);
  ~RecoTrainingNtupleMaker() {};

  void Run();

};

#endif
