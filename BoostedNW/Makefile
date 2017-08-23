ROOTLIBS  = `root-config --libs`
ROOTFLAGS = `root-config --cflags --glibs`
#SPECIALFLAGS= --exceptions
CFLAGS = -I.
RCXX =$(DELPLIBS) $(EXTLIBS) $(CFLAGS) $(ROOTFLAGS) $(ROOTLIBS) -lTMVA


DBG = -g -O2
CC = g++ $(DBG) $(RCXX)
SRC = src
UTILS= utils
INCLUDE = include
LIBS= lib
BIN= bin


all: NeutrinoWeighter BNWReconstruction TMVAClassification main reconstruct

NeutrinoWeighter: $(SRC)/NeutrinoWeighter.cxx 
	$(CC) -o $(LIBS)/NeutrinoWeighter.o -c $(SRC)/NeutrinoWeighter.cxx
	
BNWReconstruction: $(SRC)/BNWReconstruction.cxx $(INCLUDE)/NeutrinoWeighter.h 
	$(CC) -o $(LIBS)/BNWReconstruction.o -c $(SRC)/BNWReconstruction.cxx

TMVAClassification: $(SRC)/TMVAClassification.cxx
	$(CC) -o $(LIBS)/TMVAClassification.o -c $(SRC)/TMVAClassification.cxx

main: $(SRC)/main.cxx
	$(CC) -o $(BIN)/BDTTrainer $(SRC)/main.cxx $(LIBS)/TMVAClassification.o
reconstruct: $(SRC)/reconstruct.cxx
	$(CC) -o $(BIN)/BNWReco $(SRC)/reconstruct.cxx $(LIBS)/BNWReconstruction.o $(LIBS)/NeutrinoWeighter.o 
clean:
	@rm $(EXEC)
