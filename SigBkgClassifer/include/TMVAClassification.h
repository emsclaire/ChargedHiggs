#include "include/RootIncludes.h"


class TMVAClassification
{
private:

  void ReadInputFiles();
  void GetSigBkgTrees();
  void Book();
  void TrainAndTest();

  std::string m_inputSigFileName;
  std::string m_inputBkgFileName;
  std::string m_typeBDT;
  std::string m_region;
  std::string m_outputTMVA;
  std::string m_outputTrain;

  TFile *m_inputSigFile;
  TFile *m_inputBkgFile;
  TTree *m_sigTree;
  TTree *m_bkgTree;

  TFile *m_outputFile;

  TMVA::Factory *m_factory;

public:

  TMVAClassification(std::string inputSigFileName, std::string inputBkgFileName, std::string typeBDT, std::string region);
  ~TMVAClassification() {};
	void Run();

};
