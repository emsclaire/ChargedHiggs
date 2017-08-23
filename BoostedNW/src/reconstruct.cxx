#include "include/BNWReconstruction.h"
#include "include/RootIncludes.h"

int main(int argc, char* argv[]) {

  if ( argc < 4 || argc > 5 ) abort();
  else if (argc == 4) {

    const char* inputFile    = argv[1];
    const char* outputFile   = argv[2];
    const char* typeBDT      = argv[3]; //m300, m900

    BNWReconstruction* reconstruct = new BNWReconstruction(inputFile, outputFile, typeBDT);
    reconstruct->Run();
    delete reconstruct;

    std::cout << "** Exiting..." << std::endl;

  }
  return 0;
}
