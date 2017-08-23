#include "include/Reconstruction.h"
#include "include/RecoLevel.h"
#include "include/PartonLevel.h"
#include "include/AnalysisIncludes.h"

//---------------------
// convert int to string
//---------------------
std::string IntToString(int a) {
  std::ostringstream temp;
  temp << a;
  return temp.str();
}

//---------------------
// Get cross-section for given sample
//---------------------
float GetNormalisationToLumi(std::string input, float lumi, int nEvents) {
  float crossSection;

  if (input.find("m-300") != std::string::npos) {
    crossSection = 0.00169*1000.;
  }
  else if (input.find("m400.") != std::string::npos) {
    crossSection = 0.000878*1000;
  }
  else if (input.find("m500.") != std::string::npos) {
    crossSection = 0.00048039*1000;
  }
  else if (input.find("m600.") != std::string::npos) {
    crossSection = 0.00027465*1000;
  }
  else if (input.find("m700.") != std::string::npos) {
    crossSection = 0.00016344*1000;
  }
  else if (input.find("m800.") != std::string::npos) {
    crossSection = 0.00010023*1000;
  }
  else if (input.find("m-900") != std::string::npos) {
    crossSection = 0.0000631*1000.;
  }
  else if (input.find("ttbb") != std::string::npos) {
    crossSection = 0.03373*1000.;
  }
  else if (input.find("ttcc") != std::string::npos) {
    crossSection = 0.02660*1000.;
  }
  else if (input.find("ttlight") != std::string::npos) {
    crossSection = 0.12672*1000.;
  }

  std::cout << "Cross-section: " << crossSection << " fb" << std::endl;

  float evts_over_xsec = (float)nEvents/crossSection;
  return (lumi/evts_over_xsec);
}


//---------------------
// for run i, add 25000*(i-1) to event number
//---------------------
int EvtMultiplier(std::string inFileName) {
  int multiplier = -1;

  bool foundRun = false;
  int i = 0;
  while (!foundRun) {
    i++;
    std::string runNum;
    if (i >= 10) runNum = "run-"+IntToString(i);
    else runNum = "run-0"+IntToString(i);
    // determine which sample we are in and set muliplier for events
    if (inFileName.find(runNum) != std::string::npos) {
      multiplier = 25000*(i-1);
      foundRun = true;
    }
  }
  return multiplier;
}

//---------------------
// Sort objects by pt
//---------------------

template <typename T>
bool SortByPt_main(T obj1, T obj2) {
  return ( obj1->PT > obj2->PT );
}

template <typename T1, typename T2>
bool SortByPt_main(std::pair<T1,T2> obj1, std::pair<T1,T2> obj2) {
  return (obj1.second)->PT > (obj2.second)->PT;
}

bool SortLorentzByPt_main(std::pair<std::pair<std::string,int>,TLorentzVector> obj1, std::pair<std::pair<std::string,int>,TLorentzVector> obj2) {
  return ( obj1.second.Pt() > obj2.second.Pt() );
}

//---------------------
// DeltaR (taking into account that it should be between 0 and 2pi)
//---------------------
float DeltaR(float eta1, float eta2, float phi1, float phi2) {
  float result = -100000.;
  float deltaEta = eta1 - eta2;
  float deltaPhi = phi1 - phi2;
  double kPI = TMath::Pi();
  double kTWOPI = 2.*kPI;
  while (deltaPhi >= kPI) deltaPhi -= kTWOPI;
  while (deltaPhi < -kPI) deltaPhi += kTWOPI;
  result = sqrt( (deltaEta*deltaEta) + (deltaPhi*deltaPhi) );
  return result;
}

//---------------------
// get hard scatter particles and store those passing cuts
//---------------------

// find mother of particle (as long as has specified pid)
int Mother(GenParticle* particle, TClonesArray* branchParticle, int pid) {
  int result = 0;
  // should not need, in general, to go back more than 5 particles to find mother of interest
  for (int i = 0; i < 5; i++) {
    if (particle->M1+1 != 0) {
      GenParticle* mother = (GenParticle*)branchParticle->At(particle->M1);
      // use abs for charged Higgs because interested in whether it is charged Higgs, but not the sign of Hc
      if (pid == 37) {
        if (abs(mother->PID) != pid) particle = mother;
        else {
          result = particle->M1+1;
          break;
        }
      }
      else {
        if ( mother->PID != pid ) particle = mother;
        else {
          result = particle->M1+1;
          break;
        }
      }
    }
  }
  return result;
}

// swap particles if pt is lower in "1" than "2"
void SwapParticles(std::map<std::string,GenParticle*>& genPar, std::string parType) {
  if (genPar[parType+"FromTop1"]->PT < genPar[parType+"FromTop2"]->PT) {
    GenParticle* temp = genPar[parType+"FromTop1"];
    genPar[parType+"FromTop1"] = genPar[parType+"FromTop2"];
    genPar[parType+"FromTop2"] = temp;
  }
  return;
}

// get hard scatter truth particles
void GetHardScatterParticles(bool m_sig, TClonesArray* branchParticle, TIter* itGenParticle, std::map<std::string,GenParticle*>& GenParticles, std::map<std::string,TLorentzVector>& GenParticlesLorentz, std::vector<std::pair<std::string,GenParticle*> >& bQuarks, std::vector<std::pair<std::string,TLorentzVector> >& bQuarksLorentz) {

  // index denotes row in which the mother is found (I think)
  int topIndex = -1;
  int higgsIndex = -1;
  int lep_topIndex = -1;
  int lep_higgsIndex = -1;
  int nu_topIndex = -1;
  int nu_higgsIndex = -1;

  int chargeMultiplier = 0;

  int n_bFromG = 0;
  int n_bFromTop = 0;
  int n_lepFromTop = 0;
  int n_nuFromTop = 0;

  GenParticle* particle;

  while ( ( particle = (GenParticle*)itGenParticle->Next() ) ) {
    if ( particle->Status == 3 ) {

      // charged Higgs
      if (abs(particle->PID) == 37) {
        GenParticles["ChargedHiggs"] = particle;
      }

      // top
      else if (particle->PID == 6) {
        GenParticles["top"] = particle;
      }

      // tbar
      else if (particle->PID == -6) {
        GenParticles["tbar"] = particle;
      }

      // b-partons
      else if (abs(particle->PID) == 5) {
        // as negative PID is anti b-quark, need negative top PID (i.e. anti-top)
        if ( particle->PID < 0 ) chargeMultiplier = -1;
        else chargeMultiplier = +1;
        topIndex = Mother(particle, branchParticle, 6*chargeMultiplier);
        higgsIndex = Mother(particle, branchParticle, 37);

        // Hc -> t -> b (signal only)
        if (higgsIndex > 0 && topIndex > 0) {
          GenParticles["bFromTopAndHiggs"] = particle;
          GenParticlesLorentz["bFromTopAndHiggs"] = particle->P4();
        }
        // g -> t -> b
        else if (higgsIndex <= 0 && topIndex > 0) {
          // if signal, this is the only g -> t -> b
          if(m_sig) {
            GenParticles["bFromTop"] = particle;
            GenParticlesLorentz["bFromTop"] = particle->P4();
          }
          // if background, there are two g -> t -> b decays (with opposite charge)
          else {
            n_bFromTop++;
            if (n_bFromTop == 1) {
              GenParticles["bFromTop1"] = particle;
              GenParticlesLorentz["bFromTop1"] = particle->P4();
            }
            else if (n_bFromTop == 2) {
              GenParticles["bFromTop2"] = particle;
              GenParticlesLorentz["bFromTop2"] = particle->P4();
            }
          }
        }
        // Hc -> b (signal only)
        else if (higgsIndex > 0 && topIndex <= 0) {
          GenParticles["bFromHiggs"] = particle;
          GenParticlesLorentz["bFromHiggs"] = particle->P4();
        }
        // g -> b
        else {
          // if signal, this is the only g -> b
          if(m_sig) {
            GenParticles["bFromGluon"] = particle;
            GenParticlesLorentz["bFromGluon"] = particle->P4();
          }
          // in ttbb background, gluon splitting means there is g -> b bbar i.e. two
          else {
            n_bFromG++;
            if (n_bFromG == 1) {
              GenParticles["bFromGluon1"] = particle;
              GenParticlesLorentz["bFromGluon1"] = particle->P4();
            }
            else if (n_bFromG == 2) {
              GenParticles["bFromGluon2"] = particle;
              GenParticlesLorentz["bFromGluon2"] = particle->P4();
            }
          }
        }
      }

      // leptons
      else if (abs(particle->PID) == 11 || abs(particle->PID) == 13) {
        // as negative PID is a positively charged lepton (i.e. anti-lepton), need positive top PID (top not anti-top)
        if ( particle->PID < 0 ) chargeMultiplier = +1;
        else chargeMultiplier = -1;
        lep_topIndex = Mother(particle, branchParticle, 6*chargeMultiplier);
        lep_higgsIndex = Mother(particle, branchParticle, 37);

        // Hc -> t -> l (signal only)
        if (lep_higgsIndex > 0 && lep_topIndex > 0) {
          GenParticles["lepFromTopAndHiggs"] = particle;
          GenParticlesLorentz["lepFromTopAndHiggs"] = particle->P4();
        }
        // g -> t -> l
        else if ( lep_higgsIndex <= 0 && lep_topIndex > 0) {
          // if signal, there is only one g -> t -> l decay
          if(m_sig) {
            GenParticles["lepFromTop"] = particle;
            GenParticlesLorentz["lepFromTop"] = particle->P4();
          }
          // in background, have two g -> t -> l decays (with opposite charge)
          else {
            n_lepFromTop++;
            if (n_lepFromTop == 1) {
              GenParticles["lepFromTop1"] = particle;
              GenParticlesLorentz["lepFromTop1"] = particle->P4();
            }
            else if (n_lepFromTop == 2) {
              GenParticles["lepFromTop2"] = particle;
              GenParticlesLorentz["lepFromTop2"] = particle->P4();
            }
          }
        }
        else { std::cout << "Cannot have lep_higgsIndex = " << lep_higgsIndex << " and lep_topIndex = " << lep_topIndex << std::endl; }
      }

      // neutrinos
      else if (abs(particle->PID) == 12 || abs(particle->PID) == 14 ) {
        // as negative PID is antineutrino, need negative top PID (i.e. anti-top)
        if ( particle->PID < 0 ) chargeMultiplier = -1;
        else chargeMultiplier = +1;
        nu_topIndex = Mother(particle, branchParticle, 6*chargeMultiplier);
        nu_higgsIndex = Mother(particle, branchParticle, 37);

        // Hc -> t -> nu (signal only)
        if (nu_higgsIndex > 0 && nu_topIndex > 0) {
          GenParticles["nuFromTopAndHiggs"] = particle;
          GenParticlesLorentz["nuFromTopAndHiggs"] = particle->P4();
        }
        // g -> t -> nu
        else if ( nu_higgsIndex <= 0 && nu_topIndex > 0) {
          // if signal, have only one g -> t -> nu decay
          if(m_sig) {
            GenParticles["nuFromTop"] = particle;
            GenParticlesLorentz["nuFromTop"] = particle->P4();
          }
          // if background, have two g -> t -> nu decays (with opposite charge)
          else {
            n_nuFromTop++;
            if (n_nuFromTop == 1) {
              GenParticles["nuFromTop1"] = particle;
              GenParticlesLorentz["nuFromTop1"] = particle->P4();
            }
            else if (n_nuFromTop == 2) {
              GenParticles["nuFromTop2"] = particle;
              GenParticlesLorentz["nuFromTop2"] = particle->P4();
            }
          }
        }
        else { std::cout << "Cannot have nu_higgsIndex = " << nu_higgsIndex << " and nu_topIndex = " << nu_topIndex << std::endl; }
      }
    }
  }

  std::map<std::string,GenParticle*>::iterator gen_itr = GenParticles.begin();
  std::map<std::string,GenParticle*>::iterator gen_end = GenParticles.end();
  for (; gen_itr != gen_end; ++gen_itr) {
    if ((gen_itr->first).find("bFrom") != std::string::npos) {
      bQuarks.push_back(*gen_itr);
      bQuarksLorentz.push_back(std::make_pair((*gen_itr).first,(*gen_itr).second->P4()));
    }
  }

  // sort b-partons by pt regardless of origin of b-partons
  std::sort(bQuarks.begin(), bQuarks.end(),
      [](const std::pair<std::string,GenParticle*>& q1, const std::pair<std::string,GenParticle*>& q2) {
        return ( (q1.second)->PT > (q2.second)->PT );
      });
  std::sort(bQuarksLorentz.begin(), bQuarksLorentz.end(),
      [](const std::pair<std::string,TLorentzVector> qL1, const std::pair<std::string,TLorentzVector> qL2) {
        return ( qL1.second.Pt() > qL2.second.Pt() );
      });

  return;
}

void GetAllJets(TIter* itJet, std::vector<Jet*>& v_AllJets, std::vector<TLorentzVector>& v_AllJetsLorentz, int& num_btags, int& m_passJetKine, int& m_passBtag, bool& m_HasPassedJetCuts) {
  // apply cuts on jets (pt, eta)
  int num_jetkine = 0;
  num_btags = 0;
  Jet* jet;
  while ( ( jet = (Jet*)itJet->Next() ) ) {
    if( jet->PT >= 30 && fabs(jet->Eta) < 2.4 ) {
      num_jetkine++;
      v_AllJets.push_back(jet);
      v_AllJetsLorentz.push_back(jet->P4());
      if (jet->BTag == 1 ) num_btags++;
    }
  }

  if (num_jetkine >= 4) {
    m_passJetKine++;
    if (num_btags >= 2) {
      m_passBtag++;
      m_HasPassedJetCuts = true;
    }
  }

  return;
}

//---------------------
// Apply lepton cuts
//---------------------
template <typename T, typename S>
void GetDeltaRAndInvariantMass(T vec1, S vec2, float &deltaR, float &minv) {
  deltaR = DeltaR(vec1->Eta, vec2->Eta, vec1->Phi, vec2->Phi);
  TLorentzVector lorentzVec1 = vec1->P4();
  TLorentzVector lorentzVec2 = vec2->P4();
  minv = (lorentzVec1 + lorentzVec2).M();
  return;
}

//---------------------
// get reco leptons and store those passing cuts
//---------------------
void GetRecoLeptons(int evt, TIter* itElectron, TIter* itMuon, std::vector<std::pair<std::pair<std::string,int>, TLorentzVector> >& RecoLeptonsLorentz, int& m_passLepKine, int& m_passLepOppCharge, int& m_passIsoMinv, bool& m_HasPassedLepCuts) {
  //  std::cout << "GetRecoLeptons()" << std::endl;

  std::vector<std::pair<std::pair<std::string,int>, TLorentzVector> > v_allLeptons;

  // apply cuts on leptons (pt, eta)
  Electron* elec;
  while ( ( elec = (Electron*)itElectron->Next() ) ) {
    v_allLeptons.push_back(std::make_pair(std::make_pair("e",elec->Charge),elec->P4()));
  }
  Muon* muon;
  while ( ( muon = (Muon*)itMuon->Next() ) ) {
    v_allLeptons.push_back(std::make_pair(std::make_pair("m",muon->Charge),muon->P4()));
  }

  // sort by pt
  std::sort(v_allLeptons.begin(), v_allLeptons.end(), SortLorentzByPt_main);

  float dR;
  float minv;
  bool oppositeCharge = false;

  if (v_allLeptons.size() == 2) {
    if (v_allLeptons[0].second.Pt() >= 20 && fabs(v_allLeptons[0].second.Eta()) < 2.5) {
      if (v_allLeptons[1].second.Pt() >= 17 && fabs(v_allLeptons[1].second.Eta()) < 2.5) {
        m_passLepKine++;
        if (v_allLeptons[0].first.second == -(v_allLeptons[1].first.second)) {
          m_passLepOppCharge++;
          dR = (v_allLeptons[0].second).DeltaR(v_allLeptons[1].second);
          minv = (v_allLeptons[0].second + v_allLeptons[1].second).M();
          if (dR > 0.4 && minv > 12) {
            // if same sign
            // std::cout << dR << " " << minv << " " << v_allLeptons[0].first.first << " " << v_allLeptons[1].first.first << std::endl;
            if (v_allLeptons[0].first.first.compare(v_allLeptons[1].first.first) == 0) {
              if ((fabs(minv - 91.1876) > 10)) {
                m_passIsoMinv++;
                m_HasPassedLepCuts = true;
                RecoLeptonsLorentz = v_allLeptons;
              }
            }
            else {
              m_passIsoMinv++;
              m_HasPassedLepCuts = true;
              RecoLeptonsLorentz = v_allLeptons;
            }
          }
        }
      }
    }
  }

  return;
}

//---------------------
// get MET and store if passes cuts
//---------------------
MissingET* GetMET(TIter* itMissingET, int& m_passMET) {
  // apply cuts on MET
  MissingET* missET;
  while ( ( missET = (MissingET*)itMissingET->Next() ) ) {
    if( missET->MET >= 40 ) {
      m_passMET++;
      return missET;
    }
    else return NULL;
  }
}

int main(int argc, char* argv[]) {

  gROOT->SetBatch(kTRUE);
  TStyle *atlasStyle = (TStyle*)gROOT->GetListOfSpecials()->FindObject("ATLAS");
  // Take the submit directory and inputFilePath from the input if provided:
  if ( argc < 5 || argc > 6 ) abort();
  else if (argc == 5) {
    // set input values
    const char* inputFileName  = argv[1];
    std::string inFileName     = argv[1];
    const char* outputFileName = argv[2];
    bool m_sig       = (strcmp(argv[3],"sig") == 0);
    bool m_PrintInfo = (strcmp(argv[4],"print") == 0);
    int evt_multiplier = 0 ;//EvtMultiplier(inFileName);

    // the chain retrieves the tree contained in the rootfile called in argument.
    TChain* chain = new TChain("Delphes");
    chain->Add(inputFileName);

    // file for results
    TFile* ofile = new TFile(outputFileName, "recreate");
    TDirectory *saveDirectory = gDirectory;
    saveDirectory->cd();
    // treereader to run over the chain.
    ExRootTreeReader* treeReader = new ExRootTreeReader(chain);
    // branches in the tree
    TClonesArray* branchParticle  = treeReader->UseBranch("Particle");
    TClonesArray* branchGenJet    = treeReader->UseBranch("GenJet");
    TClonesArray* branchElectron  = treeReader->UseBranch("Electron");
    TClonesArray* branchMuon      = treeReader->UseBranch("Muon");
    TClonesArray* branchJet       = treeReader->UseBranch("Jet");
    TClonesArray* branchMissingET = treeReader->UseBranch("MissingET");
    // number of events in tree
    Long64_t allEntries = treeReader->GetEntries();
    // iterators over branches
    TIter* itGenParticle = new TIter(branchParticle);
    TIter* itGenJet      = new TIter(branchGenJet);
    TIter* itElectron    = new TIter(branchElectron);
    TIter* itMuon        = new TIter(branchMuon);
    TIter* itJet         = new TIter(branchJet);
    TIter* itMissingET   = new TIter(branchMissingET);

    // normalise cut flow event numbers to luminosity
    float luminosity = 100.;
    float m_normalisation = GetNormalisationToLumi(inFileName, luminosity, allEntries);

    //-------------------
    // Create objects and book histograms
    //-------------------
    Reconstruction* reco = new RecoLevel(m_normalisation, m_sig, m_PrintInfo);
    reco->BookTreeBranches();

    Reconstruction* parton = new PartonLevel(m_normalisation, m_sig, m_PrintInfo);
    parton->BookTreeBranches();

    int m_correctHardScatter = 0;
    int m_passJetKine = 0;
    int m_passBtag = 0;
    int m_passLepKine = 0;
    int m_passLepOppCharge = 0;
    int m_passIsoMinv = 0;
    int m_passMET = 0;

    //-------------------
    // loop through events
    //-------------------
    for ( int evt = 0; evt < allEntries; evt++ ) {

      std::map<std::string,GenParticle*> in_GenParticles;
      std::map<std::string,TLorentzVector> in_GenParticlesLorentz;
      std::vector<std::pair<std::string,GenParticle*> > in_bQuarks;
      std::vector<std::pair<std::string,TLorentzVector> > in_bQuarksLorentz;
      std::vector<Jet*> in_AllJets;
      std::vector<TLorentzVector> in_AllJetsLorentz;
      int numberOfBTags = 0;
      std::vector<std::pair<std::pair<std::string,int>, TLorentzVector> > in_RecoLeptonsLorentz;
      MissingET* in_MissEt = NULL;

      bool m_HasPassedJetCuts = false;
      bool m_HasPassedLepCuts = false;
      bool m_HasPassedMETCut = false;

      if (m_PrintInfo) std::cout << "Event " << evt << std::endl;
      else {
        if (evt%250 == 0) std::cout << "Event " << evt << std::endl;
      }
      treeReader->ReadEntry(evt);
      // reset iterators over branches
      itGenParticle->Reset();
      itGenJet->Reset();
      itElectron->Reset();
      itMuon->Reset();
      itJet->Reset();
      itMissingET->Reset();

      //-------------------
      // get objects at parton-level
      //-------------------
      GetHardScatterParticles(m_sig, branchParticle, itGenParticle, in_GenParticles, in_GenParticlesLorentz, in_bQuarks, in_bQuarksLorentz);

      //------------------
      // If background, order b/lep/nu from top by pt, and then b/lep/nu from gluon by pt
      // i.e. 0th element: objFromTop with highest pt
      //      1st element: objFromTop with lowest pt
      //      2nd element: objFromGluon with highest pt
      //      3rd element: objFromGluon with lowest pt
      //------------------
      if(!m_sig) {
        SwapParticles(in_GenParticles, "b");
        SwapParticles(in_GenParticles, "lep");
        SwapParticles(in_GenParticles, "nu");
      }
      //------------------
      // Insert hard scatter particles into class objects
      //------------------
      reco->SetHardScatterParticles(in_GenParticles, in_GenParticlesLorentz, in_bQuarks, in_bQuarksLorentz);
      parton->SetHardScatterParticles(in_GenParticles, in_GenParticlesLorentz, in_bQuarks, in_bQuarksLorentz);

      //------------------
      // Get and insert jets into class objects, provided >= 4 jets pass cuts
      //------------------
      if ( branchJet->GetEntriesFast() >= 4 ) {
        GetAllJets(itJet, in_AllJets, in_AllJetsLorentz, numberOfBTags, m_passJetKine, m_passBtag, m_HasPassedJetCuts);
        reco->SetBJets(in_AllJets, in_AllJetsLorentz, numberOfBTags, 2); // see Reconstruction.cxx
      }
      //------------------
      // Get and insert leptons into class objects, provided exactly 2 leptons are generated and pass cuts
      //------------------
      if (m_HasPassedJetCuts) {
        if ( (branchElectron->GetEntriesFast()+branchMuon->GetEntriesFast()) == 2) {
          GetRecoLeptons(evt, itElectron, itMuon, in_RecoLeptonsLorentz, m_passLepKine, m_passLepOppCharge, m_passIsoMinv, m_HasPassedLepCuts);
          reco->SetRecoLeptons(in_RecoLeptonsLorentz); // see Reconstruction.cxx
        }

        //------------------
        // Get and insert MET into class objects, provided pass cuts
        //------------------
        if (m_HasPassedLepCuts) {
          in_MissEt = GetMET(itMissingET, m_passMET);
          m_HasPassedMETCut = (in_MissEt != NULL);
          reco->SetMET(in_MissEt); // see Reconstruction.cxx
        }

      }
      //------------------
      // reconstruct event
      //------------------
      bool m_analyse = false;
      if (m_sig) {
        m_analyse = (in_GenParticles["bFromTop"] && in_GenParticles["bFromTopAndHiggs"] && in_GenParticles["bFromHiggs"] && in_GenParticles["bFromGluon"] && in_GenParticles["lepFromTop"] && in_GenParticles["lepFromTopAndHiggs"] && in_GenParticles["nuFromTop"] && in_GenParticles["nuFromTopAndHiggs"]);
      }
      else m_analyse = true;

      if (m_analyse) {
        reco->Event(inputFileName, outputFileName, evt+evt_multiplier, numberOfBTags); // see RecoLevel.cxx
        parton->Event(inputFileName, outputFileName, evt+evt_multiplier, -1); // see PartonLevel.cxx
      }

      reco->ClearEventVectors();
      parton->ClearEventVectors();

      //------------------
      // Clear maps and vectors
      //------------------
      in_GenParticles.clear();
      in_GenParticlesLorentz.clear();
      in_bQuarks.clear();
      in_bQuarksLorentz.clear();
      in_AllJets.clear();
      in_AllJetsLorentz.clear();
      in_RecoLeptonsLorentz.clear();

    }

    std::cout << std::endl;
    std::cout << "jet cuts:" << std::endl;
    std::cout << "   pt & eta: From " << allEntries*m_normalisation << " to " << m_passJetKine*m_normalisation << std::endl;
    std::cout << "   btag: From " << m_passJetKine*m_normalisation << " to " << m_passBtag*m_normalisation << std::endl;
    std::cout << "lep cuts:" << std::endl;
    std::cout << "   pt & eta: From " << m_passBtag*m_normalisation << " to " << m_passLepKine*m_normalisation << std::endl;
    std::cout << "   oppCharge: From " << m_passLepKine*m_normalisation << " to " << m_passLepOppCharge*m_normalisation << std::endl;
    std::cout << "   iso & minv: From " << m_passLepOppCharge*m_normalisation << " to " << m_passIsoMinv*m_normalisation << std::endl;
    std::cout << "MET cut:" << std::endl;
    std::cout << "   E: From " << m_passIsoMinv*m_normalisation << " to " << m_passMET*m_normalisation << ": " << std::endl;
    std::cout << std::endl;

    //-------------------
    // Write trees in each case
    //-------------------
    TTree* recoTree   = reco->GetTree();
    TTree* partonTree = parton->GetTree();
    recoTree->Write();
    partonTree->Write();
    reco->DeleteTrees();
    parton->DeleteTrees();

    reco->PrintInfo();
    parton->PrintInfo();

    reco->ClearTreeVectors();
    parton->ClearTreeVectors();

    //-------------------
    // Delete allocated memory
    //-------------------
    ofile->Close();
    delete ofile;
    delete chain;
    delete treeReader;
    delete itGenParticle;
    delete itGenJet;
    delete itElectron;
    delete itMuon;
    delete itJet;
    delete itMissingET;

    delete reco;
    delete parton;

    std::cout << "** Exiting..." << std::endl;

  }
  return 0;
}
