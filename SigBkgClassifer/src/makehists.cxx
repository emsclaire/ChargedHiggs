#include "include/LimitsHistMaker.h"
#include "include/RootIncludes.h"

int main(int argc, char* argv[]) {

  if (argc < 5 || argc > 6) abort();
  else if (argc == 5) {

    const char* inputFile    = argv[1];
    const char* outputFile   = argv[2];
    const char* typeBDT      = argv[3]; //m300, m900

    std::string s_bdtWeightCut = argv[4];
    bool haveBDTWeightCut = false;
    haveBDTWeightCut = !(s_bdtWeightCut.compare("NoCut") == 0);
    float bdtWeightCutValue = -99999.;
    if (haveBDTWeightCut) bdtWeightCutValue = std::stof(s_bdtWeightCut);

    float normalisation = 1.;
    std::string s_inputFile(inputFile);

    if (s_inputFile.find("m300.") != std::string::npos) {
      normalisation = 100. / (25000./(0.00169*1000));
    }
    else if (s_inputFile.find("m400.") != std::string::npos) {
      normalisation = 100. / (25000./(0.00087821*1000));
    }
    else if (s_inputFile.find("m500.") != std::string::npos) {
      normalisation = 100. / (25000./(0.00048039*1000));
    }
    else if (s_inputFile.find("m600.") != std::string::npos) {
      normalisation = 100. / (25000./(0.00027465*1000));
    }
    else if (s_inputFile.find("m700.") != std::string::npos) {
      normalisation = 100. / (25000./(0.00016344*1000));
    }
    else if (s_inputFile.find("m800.") != std::string::npos) {
      normalisation = 100. / (25000./(0.00010023*1000));
    }
    else if (s_inputFile.find("m900.") != std::string::npos) {
      normalisation = 100. / (25000./(0.0000631*1000));
    }
    else if (s_inputFile.find("ttbb") != std::string::npos) {
      normalisation = 100. / (25000./(0.03373*1000));
    }
    else if (s_inputFile.find("ttcc") != std::string::npos) {
      normalisation = 100. / (25000./(0.0266*1000));
    }
    else if (s_inputFile.find("ttlight") != std::string::npos) {
      normalisation = 100. / (25000./(0.12672*1000));
    }

    LimitsHistMaker* hmaker = new LimitsHistMaker(inputFile, outputFile, typeBDT, haveBDTWeightCut, bdtWeightCutValue, normalisation);
    hmaker->Run();
    delete hmaker;

    std::cout << "** Exiting..." << std::endl;

  }
  return 0;
}
