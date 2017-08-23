#ifndef PartonLevel_h
#define PartonLevel_h
#include "include/Reconstruction.h"
#include "include/AnalysisIncludes.h"

class PartonLevel: public Reconstruction {

private:
  // set tree variables to default value
  void InitialiseEventVariables();
  // print event information
  void PrintEventInfo();

public:
  // constructor and destructor
  PartonLevel(float in_norm, bool in_sig, bool in_print) : Reconstruction(in_norm,in_sig,in_print) {
    deltaR_bquarks.resize(6);
    v_bPartons_Lorentz.resize(4);
    v_bPartons_Origin.resize(4);
    v_bPartons_pt.resize(4);
    v_bPartons_eta.resize(4);
    v_bPartons_phi.resize(4);
    v_bPartons_E.resize(4);
    v_bPartons_M.resize(4);
  }
  ~PartonLevel() {}

  // book trees for parton level
  void BookTreeBranches();
  // analyse one event
  void Event(const char* inputFilename, const char* outputFilename, int evtNum, int numberOfBTags);
  // print information for whole sample
  void PrintInfo();

};

#endif
