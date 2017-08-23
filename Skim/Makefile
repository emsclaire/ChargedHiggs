HOMEDIR=/afs/cern.ch/work/r/rnaranjo/public/ChargedHiggs
ROOTLIBS  = `root-config --libs`
ROOTFLAGS = `root-config --cflags --glibs`
DELPLIBS=-I$(HOMEDIR)/Delphes-3.3.2 -L$(HOMEDIR)/Delphes-3.3.2 -lDelphes
EXTLIBS=-I$(HOMEDIR)/ExRootAnalysis -L$(HOMEDIR)/ExRootAnalysis -lExRootAnalysis
#SPECIALFLAGS= --exceptions
CFLAGS = -I.
RCXX =$(DELPLIBS) $(EXTLIBS) $(CFLAGS) $(ROOTFLAGS) $(ROOTLIBS)


DBG = -g -O2
CC = g++ $(DBG) $(RCXX)
SRC = src
UTILS= utils
INCLUDE = include
LIBS= lib
BIN= bin


all: Reconstruction  PartonLevel RecoLevel main treeSplit


Reconstruction: $(SRC)/Reconstruction.cxx
	$(CC) -o $(LIBS)/Reconstruction.o -c $(SRC)/Reconstruction.cxx

PartonLevel: $(SRC)/PartonLevel.cxx $(INCLUDE)/Reconstruction.h
	$(CC) -o $(LIBS)/PartonLevel.o -c $(SRC)/PartonLevel.cxx

RecoLevel: $(SRC)/PartonLevel.cxx $(INCLUDE)/Reconstruction.h
	$(CC) -o $(LIBS)/RecoLevel.o -c $(SRC)/RecoLevel.cxx

main: $(SRC)/main.cxx
	$(CC) -o $(BIN)/Analyser $(SRC)/main.cxx $(LIBS)/Reconstruction.o $(LIBS)/PartonLevel.o $(LIBS)/RecoLevel.o

treeSplit: $(UTILS)/treeSplit.cxx
	$(CC) -o $(BIN)/treeSplit $(UTILS)/treeSplit.cxx

clean:
	@rm DrawDict* $(EXEC)
