#include "include/TMVAClassification.h"
#include "include/RootIncludes.h"

int main(int argc, char* argv[]) {

  if ( argc < 3 || argc > 4 ) abort();
  else if (argc == 3) {

    const char* inputFile  = argv[1];
    const char* typeBDT = argv[2];

    TMVAClassification* trainRecoBDT = new TMVAClassification(inputFile,  typeBDT);
    trainRecoBDT->Run();
    delete trainRecoBDT;

    std::cout << "** Exiting..." << std::endl;
  }
  return 0;
}
