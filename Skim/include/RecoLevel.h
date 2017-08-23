#ifndef RecoLevel_h
#define RecoLevel_h
#include "include/Reconstruction.h"
#include "include/AnalysisIncludes.h"

class RecoLevel: public Reconstruction {

private:
  // set tree variables to default value
  void InitialiseEventVariables();
  // print event information
  void PrintEventInfo();

public:

  // constructor and destructor
  RecoLevel(float in_norm, bool in_sig, bool in_print) : Reconstruction(in_norm,in_sig,in_print) {

    v_bPartons_Lorentz.resize(4);
    v_bPartons_Origin.resize(4);
    v_bPartons_pt.resize(4);
    v_bPartons_eta.resize(4);
    v_bPartons_phi.resize(4);
    v_bPartons_E.resize(4);
    v_bPartons_M.resize(4);

    v_SelectedJets_Lorentz.resize(4);
    v_SelectedJets_pt.resize(4);
    v_SelectedJets_eta.resize(4);
    v_SelectedJets_phi.resize(4);
    v_SelectedJets_E.resize(4);
    v_SelectedJets_M.resize(4);
    v_SelectedJets_bTag.resize(4);

    v_bIndexMatchedToJet.resize(4);
    v_JetMatchedOrigin.resize(4);

    n_EventsPassingCuts = 0;
    n_EventsWithCorrectBQuarks = 0;
    n_EventsWithSpecificMatches = 0;
    n_EventsWith3Matches = 0;
    n_EventsWith3Matches_UsingAllJets = 0;
    n_EventsWith4Matches = 0;
    n_EventsWith4Matches_UsingAllJets = 0;
    n_Events3OrderedMatches = 0;
    n_Events4OrderedMatches = 0;
  }
  ~RecoLevel() {}

  // book trees for reco level
  void BookTreeBranches();
  // analyse one event
  void Event(const char* inputFilename, const char* outputFilename, int evtNum, int numberOfBTags);
  // print information for whole sample
  void PrintInfo();

};

#endif
