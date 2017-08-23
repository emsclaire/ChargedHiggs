ROOTLIBS  = `root-config --libs`
ROOTFLAGS = `root-config --cflags --glibs`
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


all: RecoTrainingNtupleMaker main


RecoTrainingNtupleMaker: $(SRC)/RecoTrainingNtupleMaker.cxx
	$(CC) -o $(LIBS)/RecoTrainingNtupleMaker.o -c $(SRC)/RecoTrainingNtupleMaker.cxx

main: $(SRC)/main.cxx
	$(CC) -o $(BIN)/NtupleMaker $(SRC)/main.cxx $(LIBS)/RecoTrainingNtupleMaker.o

clean:
	@rm $(EXEC)
