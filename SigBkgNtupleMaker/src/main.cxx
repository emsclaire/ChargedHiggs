#include "include/SigBkgNtupleMaker.h"
#include "include/RootIncludes.h"

int main(int argc, char* argv[]) {

  if ( argc < 3 || argc > 4 ) abort();
  else if (argc == 3) {

    std::string inputFile  = argv[1];
    std::string outputFile = argv[2];
    std::cout << "** Writing ntuple " << outputFile << "..." << std::endl;

    SigBkgNtupleMaker* ntupleMaker = new SigBkgNtupleMaker(inputFile, outputFile);
    ntupleMaker->Run();
    delete ntupleMaker;

    std::cout << "** Exiting..." << std::endl;

  }
  return 0;
}
