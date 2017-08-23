#ifndef Reconstruction_h
#define Reconstruction_h
#include "include/AnalysisIncludes.h"

// ExRootAnalysis includes
#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
// Delphes classes
#include "classes/DelphesClasses.h"


class Reconstruction {

protected:

  //----------------------
  // protected functions
  //----------------------
  std::string IntToString(int a);
  float DeltaR(float eta1, float eta2, float phi1, float phi2);
  // match jets to truth-level quarks
  void MatchJetsToQuarks();
  // alternative matching of jets to truth-level quarks
  void MatchJetsToQuarks_UsingAllJetsInEvent();
  // get positive and negative truth-level leptons
  void GetPartonPositiveNegativeLepton(TLorentzVector& leppos, TLorentzVector& lepneg);
  // get nu and nubar truth-level
  void GetPartonNeutrinoAntiNeutrino(TLorentzVector& nu, TLorentzVector& nubar);
  // set truth variables for branches
  void SetTruthVariables();
  // calculate reconstructed top/tbar invariant mass
  float CalculateTopMass(TLorentzVector b, TLorentzVector lep, TLorentzVector nu);
  // calculate reconstructed charged Higgs invariant mass
  float CalculateChargedHiggsMass(TLorentzVector b1, TLorentzVector b2, TLorentzVector lep, TLorentzVector nu);

  // set tree variables to default value
  virtual void InitialiseEventVariables() = 0;
  // print information for each event
  virtual void PrintEventInfo() = 0;

  //----------------------
  // protected variables
  //----------------------

  TTree* tree;

  // xsec of sample
  float m_norm;
  // signal or background sample
  bool m_sig;
  std::string s_level;
  // number of jets in event
  int n_totalJets;
  // number of b-tagged jets in event
  int n_btags; // up to 4
  // print out information
  bool m_print;

  // vectors to store objects
  std::map<std::string,GenParticle*> v_GenParticles;
  std::vector<Jet*> v_AllJets;
  std::vector<std::pair<int,Jet*> > v_SelectedJets;
  MissingET* m_GenMissingET;
  MissingET* m_MissingET;
  // vectors to store TLorentzVectors
  std::map<std::string,TLorentzVector> v_GenParticlesLorentz;
  std::vector<std::pair<std::pair<std::string,int>, TLorentzVector> > v_RecoLeptonsLorentz;

  std::vector<TLorentzVector> v_AllJetsLorentz;
  std::vector<std::pair<int,TLorentzVector> > v_SelectedJetsLorentz;
  // vectors for quark-jet matching
  std::vector<std::pair<std::string,GenParticle*> > v_bQuarks;
  std::vector<std::pair<std::string,TLorentzVector> > v_bQuarksLorentz;

  int eventNum;
  int passSelection;
  int numJets;
  int numBTags;
  int numMatched;
  int n_EventsPassingCuts;
  int n_EventsWithCorrectBQuarks;
  int n_EventsWithSpecificMatches;
  int n_EventsWith3Matches;
  int n_EventsWith3Matches_UsingAllJets;
  int n_EventsWith4Matches;
  int n_EventsWith4Matches_UsingAllJets;
  int n_Events3OrderedMatches;
  int n_Events4OrderedMatches;

  //------------
  // Reco variables
  //------------
  std::vector<float> v_AllJets_pt;
  std::vector<float> v_AllJets_eta;
  std::vector<float> v_AllJets_phi;
  std::vector<float> v_AllJets_E;
  std::vector<float> v_AllJets_M;
  std::vector<int> v_AllJets_bTag;

  std::vector<TLorentzVector> v_SelectedJets_Lorentz;
  std::vector<float> v_SelectedJets_pt;
  std::vector<float> v_SelectedJets_eta;
  std::vector<float> v_SelectedJets_phi;
  std::vector<float> v_SelectedJets_E;
  std::vector<float> v_SelectedJets_M;
  std::vector<int> v_SelectedJets_bTag;
  TLorentzVector delphes_poslep, delphes_neglep;
  float delphes_poslep_pt, delphes_poslep_eta, delphes_poslep_phi, delphes_poslep_E, delphes_poslep_M;
  float delphes_neglep_pt, delphes_neglep_eta, delphes_neglep_phi, delphes_neglep_E, delphes_neglep_M;
  TLorentzVector delphes_etmiss;
  float delphes_etmiss_met, delphes_etmiss_eta, delphes_etmiss_phi;

  std::vector<int> v_bIndexMatchedToJet;
  std::vector<std::pair<int,int> > v_JetQuarkPair;
  int hasRequiredSIGMatches;
  int hasRequiredBKGMatches;
  std::vector<std::string> v_JetMatchedOrigin;

  //------------
  // Truth variables
  //------------
  // charged Higgs variables
  TLorentzVector truth_Hc;
  int truth_Hc_charge;
  float truth_Hc_pt, truth_Hc_eta, truth_Hc_phi, truth_Hc_E, truth_Hc_M;
  // top variables
  TLorentzVector truth_top, truth_tbar;
  float truth_top_pt, truth_top_eta, truth_top_phi, truth_top_E, truth_top_M;
  float truth_tbar_pt, truth_tbar_eta, truth_tbar_phi, truth_tbar_E, truth_tbar_M;
  // charged lepton variables
  TLorentzVector truth_poslep, truth_neglep;
  float truth_poslep_pt, truth_poslep_eta, truth_poslep_phi, truth_poslep_E, truth_poslep_M;
  float truth_neglep_pt, truth_neglep_eta, truth_neglep_phi, truth_neglep_E, truth_neglep_M;
  // neutrino variables
  TLorentzVector truth_nu, truth_nubar;
  float truth_nu_pt, truth_nu_eta, truth_nu_phi, truth_nu_E, truth_nu_M;
  float truth_nubar_pt, truth_nubar_eta, truth_nubar_phi, truth_nubar_E, truth_nubar_M;
  // b variables for signal
  TLorentzVector bFromHiggs, bFromTopHiggs, bFromTop, bFromGluon;
  int bFromHiggs_index, bFromTopHiggs_index, bFromTop_index, bFromGluon_index;
  float bFromHiggs_pt, bFromHiggs_eta, bFromHiggs_phi, bFromHiggs_E, bFromHiggs_M;
  float bFromTopHiggs_pt, bFromTopHiggs_eta, bFromTopHiggs_phi, bFromTopHiggs_E, bFromTopHiggs_M;
  float bFromTop_pt, bFromTop_eta, bFromTop_phi, bFromTop_E, bFromTop_M;
  float bFromGluon_pt, bFromGluon_eta, bFromGluon_phi, bFromGluon_E, bFromGluon_M;
  // b variables for background
  TLorentzVector bFromTop1, bFromTop2, bFromGluon1, bFromGluon2;
  int bFromTop1_index, bFromTop2_index, bFromGluon1_index, bFromGluon2_index;
  float bFromTop1_pt, bFromTop1_eta, bFromTop1_phi, bFromTop1_E, bFromTop1_M;
  float bFromTop2_pt, bFromTop2_eta, bFromTop2_phi, bFromTop2_E, bFromTop2_M;
  float bFromGluon1_pt, bFromGluon1_eta, bFromGluon1_phi, bFromGluon1_E, bFromGluon1_M;
  float bFromGluon2_pt, bFromGluon2_eta, bFromGluon2_phi, bFromGluon2_E, bFromGluon2_M;
  // b-partons ordered by pt
  std::vector<TLorentzVector> v_bPartons_Lorentz;
  std::vector<float> v_bPartons_pt;
  std::vector<float> v_bPartons_eta;
  std::vector<float> v_bPartons_phi;
  std::vector<float> v_bPartons_E;
  std::vector<float> v_bPartons_M;
  std::vector<std::string> v_bPartons_Origin;
  // DeltaR between b-partons
  // signal
  float deltaR_bH_bt;
  float deltaR_bH_bHt;
  float deltaR_bH_bg;
  float deltaR_bHt_bt;
  float deltaR_bHt_bg;
  float deltaR_bt_bg;
  // background
  float deltaR_bt1_bt2;
  float deltaR_bt1_bg1;
  float deltaR_bt1_bg2;
  float deltaR_bt2_bg1;
  float deltaR_bt2_bg2;
  float deltaR_bg1_bg2;
  // Identifies which b-parton pairs are closest, second closest etc.
  std::string closestPair;
  std::string secondPair;
  std::string thirdPair;
  std::string fourthPair;
  std::string fifthPair;
  std::string sixthPair;
  // identifier and DeltaR pairs inserted into std::vector
  std::vector<std::pair<std::string,float> > deltaR_bquarks;

public:

  // constructor and destructor
  Reconstruction(float in_norm, bool in_sig, bool in_print): m_norm(in_norm), m_sig(in_sig), m_print(in_print) {}
  virtual ~Reconstruction() {};

  // get tree
  TTree* GetTree();
  // memory clean-up
  void DeleteTrees();
  void ClearEventVectors();
  void ClearTreeVectors();

  // get objects
  void SetHardScatterParticles(std::map<std::string,GenParticle*> in_GenParticles, std::map<std::string,TLorentzVector> in_GenParticlesLorentz, std::vector<std::pair<std::string,GenParticle*> > in_bQuarks, std::vector<std::pair<std::string,TLorentzVector> > in_bQuarksLorentz);
  // select on reco leptons
  void SetRecoLeptons(std::vector<std::pair<std::pair<std::string,int>, TLorentzVector> > in_RecoLeptonsLorentz_new);
  // select on MET
  void SetMET(MissingET* in_MissingET);
  // select on reco jets
  void SetBJets(std::vector<Jet*> in_AllJets, std::vector<TLorentzVector> in_AllJetsLorentz, int num_btags, int num_maxbtags);

  //----------------------
  // PURE VIRTUAL FUNCTIONS
  //----------------------
  // book and write histograms
  virtual void BookTreeBranches() = 0;
  // loop through events
  virtual void Event(const char* inputFilename, const char* outputFilename, int evtNum, int numberOfBTags) = 0;
  // print information for whole sample
  virtual void PrintInfo() = 0;
};

#endif
