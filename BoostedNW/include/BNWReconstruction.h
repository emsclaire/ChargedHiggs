#include "RootIncludes.h"
#include "include/NeutrinoWeighter.h"

class BNWReconstruction
{
private:
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
	void SetTMVAReader();
	void GetBDTWeight();
	void ResetBDTVariables();
	void InitialiseEventVariables();
	void BookNewTree();

	std::string m_inputFileName;

	std::string m_WeightsDir;
	std::string m_WeigthsPrefix;
	std::string m_WeightsFile;
	std::string m_typeBDT;

	std::string m_outputFileName;
	std::vector<std::string> v_inputFileList;

	TTree* m_inputTree;
	TTree* m_outputTree;

	TFile* m_inputFile;
	TFile* m_outputFile;

	TMVA::Reader* m_reader;
	




  float maxNeutrinoWeight; // parton (no origin) or reco level
  // tops
  TLorentzVector maxNW_top, maxNW_tbar;
  float maxNW_top_pt, maxNW_top_eta, maxNW_top_phi, maxNW_top_E, maxNW_top_M;
  float maxNW_tbar_pt, maxNW_tbar_eta, maxNW_tbar_phi, maxNW_tbar_E, maxNW_tbar_M;
  // bjets
  TLorentzVector maxNW_b, maxNW_bbar;
  int maxNW_b_index, maxNW_bbar_index, maxNW_bH_index, maxNW_bg_index,maxNW_btH_index, maxNW_bt_index;
  float maxNW_b_pt, maxNW_b_eta, maxNW_b_phi, maxNW_b_E, maxNW_b_M;
  float maxNW_bbar_pt, maxNW_bbar_eta, maxNW_bbar_phi, maxNW_bbar_E, maxNW_bbar_M;
  // neutrinos
  TLorentzVector maxNW_nu, maxNW_nubar;
  float maxNW_nu_pt, maxNW_nu_eta, maxNW_nu_phi, maxNW_nu_E, maxNW_nu_M;
  float maxNW_nubar_pt, maxNW_nubar_eta, maxNW_nubar_phi, maxNW_nubar_E, maxNW_nubar_M;


  int m_Reco_HcJet_index , m_Reco_tHcJet_index, m_Reco_tJet_index, m_RecogJet_index;
  float m_maxWeight, MaxBDTWeight; 
  int isReconstructed;


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
	BNWReconstruction(std::string inputFileName, std::string outputFileName, std::string typeBDT);
	~BNWReconstruction(){};
	void Run();

	
};