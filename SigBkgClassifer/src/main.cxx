#include "include/TMVAClassification.h"
#include "include/RootIncludes.h"

int main(int argc, char* argv[]) {

  if ( argc < 5 || argc > 6 ) abort();
  else if (argc == 5) {

    const char* inputSigFile  = argv[1];
    const char* inputBkgFile  = argv[2];
    const char* typeBDT = argv[3];
    const char* region = argv[4];

    TMVAClassification* ntupleMaker = new TMVAClassification(inputSigFile, inputBkgFile, typeBDT, region);
    ntupleMaker->Run();
    delete ntupleMaker;

    std::cout << "** Exiting..." << std::endl;

  }
  return 0;
}
