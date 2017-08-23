#include "include/RootIncludes.h"


class TMVAClassification
{
private:

  void ReadInputFiles();
  void SplitTrees();
  void Book();
  void TrainAndTest();
 
  std::string m_inputFileName;
  std::string m_typeBDT;
  std::string m_outputTMVA;
  std::string m_outputTrain;


  TFile *m_inputFile;
  TTree *m_bkgTree;
  TTree *m_signalTree;


  TFile *m_outputFile;

  TMVA::Factory *m_factory;



public:

	TMVAClassification(std::string inputFileName, std::string typeBDT);
  	~TMVAClassification() {};
	void Run();
	
};