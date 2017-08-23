#include <iostream>
#include <vector>
#include <algorithm>
#include <map>

#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TSystem.h>
// ExRootAnalysis includes
#include "ExRootAnalysis/ExRootClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"
// Delphes classes
#include "classes/DelphesClasses.h"


using namespace std;

int debug=0;

//forward declarations
bool IsUnique(int run, int event);


int treeSplit(TString InFileName, 
	      int     evPerFile=1000,
	      TString outFileNameBase="ntp_noDupl")
{
  gSystem->Load("libDelphes");

  // Get the ev chain.
  TFile *root_file = new TFile(InFileName);
  TTree *tree;
  TIter nextkey( root_file->GetListOfKeys() );
   TKey *key, *oldkey=0;
   while ( (key = (TKey*)nextkey())) {

      TObject *obj = key->ReadObj();
      if ( obj->IsA()->InheritsFrom( TTree::Class() ) ) {
         tree = (TTree*)obj;
         break;
      }
   }



  int evOutFile = 0;
  int fileID    = 1;

  TString outFileName = outFileNameBase+Form("_%i",fileID)+".root";
  TFile *outFile = new TFile( outFileName , "recreate");

  // New TTree for ev tree in new file
  TTree* newTree = (TTree*) tree->CloneTree(0);
  newTree->SetDirectory(outFile);
  newTree->SetMaxTreeSize((Long64_t)100000000000LL); //100GB limit

  //----------------------------------------------------
  // Loop over the ev chain and copy it to the new file
  //----------------------------------------------------
  int const nEntries = (int) tree->GetEntries();
  for (int iev=0; iev < nEntries; ++iev) {

    //if (iev % (nEntries / 100) == 0) {
    if ( iev % 10000 == 0) {
      cout << "Processing event: " 
	   << iev << " / " << nEntries << endl;
    }
    
    tree->GetEntry(iev);
    

    newTree->Fill();
    evOutFile++;

    if (evOutFile >= evPerFile) {

	//write current tree into current file
	newTree->Write();
	outFile->Close();
	outFile->Delete();
	delete outFile;
	newTree=0;
	outFile=0;

	//make a new tree, open a new file and reset counter
	evOutFile = 0;
	fileID++;

        outFileName = outFileNameBase+Form("_%i",fileID)+".root";
	outFile = new TFile( outFileName.Data() , "update");
	if (!outFile->IsOpen()) {
	  cerr << "ERROR: Cannot open file for writing.  "
	       << "  File may already exist.\n"
	       << "  This is to protect you from overwriting files.\n"
	       << "  OutFileName: " << outFileName << endl;
	  return 1;
	}
	
	// New TTree for ev tree in new file
	newTree = (TTree*) tree->CloneTree(0);
	newTree->SetDirectory(outFile);
	newTree->SetMaxTreeSize((Long64_t)100000000000LL); //100GB limit
	
     }
  } //loop over the events

  newTree->Write();
  outFile->Close();

  return 0;
}


//-----------------------------------------------------------------------------
// IsUnique
//-----------------------------------------------------------------------------
bool IsUnique(int run, int event) {

  static map< int, map<int,bool> >  themap;
  // the following line requires some explanation
  //  - themap[runno] inserts an empty map<int,bool> 
  //    into themap if it doesn't already exist
  //  - insert returns a pair the second element of which is true it
  //    the insert was sucessful and false if another entry in the way
  return themap[run].insert(pair<int,bool>(event,true)).second;
}


int help()
{
  cout << "USAGE: treeSplit <inputFileName> [<evPerFile>] [<outputFileName>]" << endl;
  cout << "\t <inputFileName> - input file with tree in it" << endl;
  cout << "\t <evPerFile> - number of events per tree in one output splitted file, default: 10000" << endl;
  cout << "\t <outputFileName> - output file name, default: ntp_noDupl.root" << endl;

  return 0;
}

//-----------------------------------------------------------------------------
// main
//-----------------------------------------------------------------------------
int main(int argc, char* argv[]) 
{
  //trivial cases --> print help
  if (argc == 1) { help(); return 0; }
  if ( argc > 1 ) {
    string arg1 = argv[1];
    if (  (arg1.compare("-h") == 0) || 
	  (arg1.compare("-H") == 0)    )  {
      help(); 
      return 0;
    }
  }

  //there is already at least one parameter now
  TString inFileName      = argv[1];
  int     evPerFile       = 1000;
  TString outFileNameBase = "ntp_noDupl";
  if (argc > 2) evPerFile       = atoi(argv[2]);
  if (argc > 3) outFileNameBase = argv[3];
  
  int rc = 0;
  
  treeSplit(inFileName, evPerFile, outFileNameBase);
  return rc;
}
