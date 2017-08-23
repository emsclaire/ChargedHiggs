#include "include/Reconstruction.h"

//---------------------
// Sort objects by pt
//---------------------

template <typename T>
bool SortByPt(T obj1, T obj2) {
  return ( obj1->PT > obj2->PT );
}

template <typename T1, typename T2>
bool SortByPt2(std::pair<T1,T2> obj1, std::pair<T1,T2> obj2) {
  return (obj1.second)->PT > (obj2.second)->PT;
}

bool SortLorentzByPt(std::pair<int,TLorentzVector> obj1, std::pair<int,TLorentzVector> obj2) {
  return ( obj1.second.Pt() > obj2.second.Pt() );
}

//---------------------
// convert int to string
//---------------------
std::string Reconstruction::IntToString(int a) {
  std::ostringstream temp;
  temp << a;
  return temp.str();
}

//---------------------
// DeltaR (taking into account that it should be between 0 and 2pi)
//---------------------
float Reconstruction::DeltaR(float eta1, float eta2, float phi1, float phi2) {
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
// PRIVATE match jet i to bquark j
//---------------------
void Reconstruction::MatchJetsToQuarks() {
  //  std::cout << "MatchJetsToQuarks()" << std::endl;

  std::vector<std::pair<int,TLorentzVector> >::iterator j_itr = v_SelectedJetsLorentz.begin();
  std::vector<std::pair<int,TLorentzVector> >::iterator j_end = v_SelectedJetsLorentz.end();

  std::vector<double> minDR;

  std::multimap<int,std::pair<int,float> > all_b_matches;
  for ( ; j_itr != j_end; ++j_itr ) {
    float DeltaRMin = 0.4;
    int i_j = std::distance(v_SelectedJetsLorentz.begin(),j_itr);
    // for each jet, iterate through the b-quarks
    std::vector<std::pair<std::string,TLorentzVector> >::iterator b_itr = v_bQuarksLorentz.begin();
    std::vector<std::pair<std::string,TLorentzVector> >::iterator b_end = v_bQuarksLorentz.end();
    for ( ; b_itr != b_end; ++b_itr ) {
      int i_b = std::distance(v_bQuarksLorentz.begin(),b_itr);
      double deltaR = (j_itr->second).DeltaR(b_itr->second);
      if (deltaR < DeltaRMin) {
        DeltaRMin = deltaR;
        all_b_matches.insert(std::pair<int,std::pair<int,float> >(i_b,std::make_pair(i_j,DeltaRMin)));
      }
    }
  }

  // find minimum DeltaR to b-quark for the given jet
  for (int i_quark = 0; i_quark < v_bQuarks.size(); i_quark++) {
    float minDR_for_b  = 0.4;
    int jetIndex_min   = 999999999;
    int quarkIndex_min = 999999999;
    std::multimap<int,std::pair<int,float> >::iterator it;
    for (it = all_b_matches.equal_range(i_quark).first; it != all_b_matches.equal_range(i_quark).second; ++it) {
      if ((it->second).second < minDR_for_b) {
        minDR_for_b    = (it->second).second;
        quarkIndex_min = i_quark;
        jetIndex_min   = (it->second.first);
      }
    }
    if (jetIndex_min != 999999999 && quarkIndex_min != 999999999) {
      v_JetQuarkPair.push_back(std::make_pair<int,int>(int(jetIndex_min),int(quarkIndex_min)));
      minDR.push_back(minDR_for_b);
    }
  }

  //---------------
  // Require specific parton-level b-quarks to have a match
  //---------------
  if (m_sig) {
    // want to have a unique match to AT LEAST bFromHiggs, bFromTopHiggs and bFromTop (i.e. bFromGluon less important)
    int nMatches_bFromHiggs = 0;
    int nMatches_bFromTopHiggs = 0;
    int nMatches_bFromTop = 0;
    std::vector<std::pair<int,int> >::iterator jqp_itr = v_JetQuarkPair.begin();
    std::vector<std::pair<int,int> >::iterator jqp_end = v_JetQuarkPair.end();
    for ( ; jqp_itr != jqp_end; ++jqp_itr) {
      if ((v_bQuarksLorentz[jqp_itr->second].first).compare("bFromHiggs") == 0) nMatches_bFromHiggs++;
      if ((v_bQuarksLorentz[jqp_itr->second].first).compare("bFromTopAndHiggs") == 0) nMatches_bFromTopHiggs++;
      if ((v_bQuarksLorentz[jqp_itr->second].first).compare("bFromTop") == 0) nMatches_bFromTop++;
    }
    if (nMatches_bFromHiggs == 1 && nMatches_bFromTopHiggs == 1 && nMatches_bFromTop == 1) hasRequiredSIGMatches = 1;
  }

  return;
}

//---------------------
// PRIVATE Obtain positive and negative parton-level leptons
//---------------------
void Reconstruction::GetPartonPositiveNegativeLepton(TLorentzVector& leppos, TLorentzVector& lepneg) {
  if (m_sig) {
    if (v_GenParticles["lepFromTop"]->Charge > 0 && v_GenParticles["lepFromTopAndHiggs"]->Charge < 0) {
      leppos = v_GenParticlesLorentz["lepFromTop"];
      lepneg = v_GenParticlesLorentz["lepFromTopAndHiggs"];
    }
    else if (v_GenParticles["lepFromTop"]->Charge < 0 && v_GenParticles["lepFromTopAndHiggs"]->Charge > 0) {
      lepneg = v_GenParticlesLorentz["lepFromTop"];
      leppos = v_GenParticlesLorentz["lepFromTopAndHiggs"];
    }
    else if (v_GenParticles["lepFromTop"]->Charge > 0 && v_GenParticles["lepFromTopAndHiggs"]->Charge > 0) {
      leppos.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      lepneg.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }

    else if (v_GenParticles["lepFromTop"]->Charge < 0 && v_GenParticles["lepFromTopAndHiggs"]->Charge < 0) {
      leppos.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      lepneg.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }
  }
  else {
    if (v_GenParticles["lepFromTop1"]->Charge > 0 && v_GenParticles["lepFromTop2"]->Charge < 0) {
      leppos = v_GenParticlesLorentz["lepFromTop1"];
      lepneg = v_GenParticlesLorentz["lepFromTop2"];
    }
    else if (v_GenParticles["lepFromTop1"]->Charge < 0 && v_GenParticles["lepFromTop2"]->Charge > 0) {
      lepneg = v_GenParticlesLorentz["lepFromTop1"];
      leppos = v_GenParticlesLorentz["lepFromTop2"];
    }
    else if (v_GenParticles["lepFromTop1"]->Charge > 0 && v_GenParticles["lepFromTop2"]->Charge > 0) {
      leppos.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      lepneg.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }

    else if (v_GenParticles["lepFromTop1"]->Charge < 0 && v_GenParticles["lepFromTop2"]->Charge < 0) {
      leppos.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      lepneg.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }
  }
  return;
}

//---------------------
// PRIVATE Obtain nu and nubar
//---------------------
void Reconstruction::GetPartonNeutrinoAntiNeutrino(TLorentzVector& nu, TLorentzVector& nubar) {
  if (m_sig) {
    if (v_GenParticles["nuFromTop"]->PID > 0 && v_GenParticles["nuFromTopAndHiggs"]->PID < 0) {
      nu    = v_GenParticlesLorentz["nuFromTop"];
      nubar = v_GenParticlesLorentz["nuFromTopAndHiggs"];
    }
    else if (v_GenParticles["nuFromTop"]->PID < 0 && v_GenParticles["nuFromTopAndHiggs"]->PID > 0) {
      nubar = v_GenParticlesLorentz["nuFromTop"];
      nu    = v_GenParticlesLorentz["nuFromTopAndHiggs"];
    }
    else if (v_GenParticles["nuFromTop"]->Charge > 0 && v_GenParticles["nuFromTopAndHiggs"]->Charge > 0) {
      nu.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      nubar.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }

    else if (v_GenParticles["nuFromTop"]->Charge < 0 && v_GenParticles["nuFromTopAndHiggs"]->Charge < 0) {
      nu.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      nubar.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }
  }
  else {
    if (v_GenParticles["nuFromTop1"]->PID > 0 && v_GenParticles["nuFromTop2"]->PID < 0) {
      nu    = v_GenParticlesLorentz["lepFromTop1"];
      nubar = v_GenParticlesLorentz["lepFromTop2"];
    }
    else if (v_GenParticles["nuFromTop1"]->Charge < 0 && v_GenParticles["nuFromTop2"]->Charge > 0) {
      nubar = v_GenParticlesLorentz["lepFromTop1"];
      nu    = v_GenParticlesLorentz["lepFromTop2"];
    }
    else if (v_GenParticles["nuFromTop1"]->Charge > 0 && v_GenParticles["nuFromTop2"]->Charge > 0) {
      nu.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      nubar.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }

    else if (v_GenParticles["nuFromTop1"]->Charge < 0 && v_GenParticles["nuFromTop2"]->Charge < 0) {
      nu.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
      nubar.SetPtEtaPhiE(-999.,-999.,-999.,-999.);
    }
  }
  return;
}

//---------------------
// PRIVATE Set truth variables for branches
//---------------------
void Reconstruction::SetTruthVariables() {

  // b-partons ordered by pt
  std::vector<std::pair<std::string,TLorentzVector> >::iterator bQuarks_itr = v_bQuarksLorentz.begin();
  std::vector<std::pair<std::string,TLorentzVector> >::iterator bQuarks_end = v_bQuarksLorentz.end();
  for ( ; bQuarks_itr != bQuarks_end; ++bQuarks_itr ) {
    int i_b = std::distance(v_bQuarksLorentz.begin(), bQuarks_itr);
    v_bPartons_Origin.at(i_b)  = bQuarks_itr->first;
    v_bPartons_Lorentz.at(i_b) = bQuarks_itr->second;
    v_bPartons_pt.at(i_b)      = v_bPartons_Lorentz.at(i_b).Pt();
    v_bPartons_eta.at(i_b)     = v_bPartons_Lorentz.at(i_b).Eta();
    v_bPartons_phi.at(i_b)     = v_bPartons_Lorentz.at(i_b).Phi();
    v_bPartons_E.at(i_b)       = v_bPartons_Lorentz.at(i_b).E();
    v_bPartons_M.at(i_b)       = v_bPartons_Lorentz.at(i_b).M();
  }

  if (m_sig) {
    // truth charged Higgs
    truth_Hc        = v_GenParticles["ChargedHiggs"]->P4();
    truth_Hc_charge = v_GenParticles["ChargedHiggs"]->Charge;
    truth_Hc_pt     = truth_Hc.Pt();
    truth_Hc_eta    = truth_Hc.Eta();
    truth_Hc_phi    = truth_Hc.Phi();
    truth_Hc_E      = truth_Hc.E();
    truth_Hc_M      = truth_Hc.M();
  }

  // truth top
  truth_top     = v_GenParticles["top"]->P4();
  truth_top_pt  = truth_top.Pt();
  truth_top_eta = truth_top.Eta();
  truth_top_phi = truth_top.Phi();
  truth_top_E   = truth_top.E();
  truth_top_M   = truth_top.M();

  // truth tbar
  truth_tbar     = v_GenParticles["tbar"]->P4();
  truth_tbar_pt  = truth_tbar.Pt();
  truth_tbar_eta = truth_tbar.Eta();
  truth_tbar_phi = truth_tbar.Phi();
  truth_tbar_E   = truth_tbar.E();
  truth_tbar_M   = truth_tbar.M();

  // Get positive and negative leptons
  GetPartonPositiveNegativeLepton(truth_poslep, truth_neglep);
  // +ve lepton
  truth_poslep_pt  = truth_poslep.Pt();
  truth_poslep_eta = truth_poslep.Eta();
  truth_poslep_phi = truth_poslep.Phi();
  truth_poslep_E   = truth_poslep.E();
  truth_poslep_M   = truth_poslep.M();
  // -ve lepton
  truth_neglep_pt  = truth_neglep.Pt();
  truth_neglep_eta = truth_neglep.Eta();
  truth_neglep_phi = truth_neglep.Phi();
  truth_neglep_E   = truth_neglep.E();
  truth_neglep_M   = truth_neglep.M();

  // Get nu and nubar
  GetPartonNeutrinoAntiNeutrino(truth_nu, truth_nubar);
  // truth nu
  truth_nu_pt  = truth_nu.Pt();
  truth_nu_eta = truth_nu.Eta();
  truth_nu_phi = truth_nu.Phi();
  truth_nu_E   = truth_nu.E();
  truth_nu_M   = truth_nu.M();
  // truth nubar
  truth_nubar_pt  = truth_nubar.Pt();
  truth_nubar_eta = truth_nubar.Eta();
  truth_nubar_phi = truth_nubar.Phi();
  truth_nubar_E   = truth_nubar.E();
  truth_nubar_M   = truth_nubar.M();

  return;
}

//---------------------
// PRIVATE Calculate top/tbar mass
//---------------------
float Reconstruction::CalculateTopMass(TLorentzVector b, TLorentzVector lep, TLorentzVector nu) {
  TLorentzVector sum = b+lep+nu;
  return sum.M();
}

//---------------------
// PRIVATE Calculate charged Higgs mass
//---------------------
float Reconstruction::CalculateChargedHiggsMass(TLorentzVector b1, TLorentzVector b2, TLorentzVector lep, TLorentzVector nu) {
  TLorentzVector sum = b1+lep+nu+b2;
  return sum.M();
}

//---------------------
// PUBLIC Clear vectors for the next event
//---------------------
void Reconstruction::ClearEventVectors() {

  // objects
  v_GenParticles.clear();
  v_GenParticlesLorentz.clear();
  v_RecoLeptonsLorentz.clear();
  v_AllJets.clear();
  v_AllJetsLorentz.clear();
  v_AllJets_pt.clear();
  v_AllJets_eta.clear();
  v_AllJets_phi.clear();
  v_AllJets_E.clear();
  v_AllJets_M.clear();
  v_AllJets_bTag.clear();
  v_SelectedJets.clear();
  v_SelectedJetsLorentz.clear();
  v_bQuarks.clear();
  v_bQuarksLorentz.clear();

  v_JetQuarkPair.clear();

  return;
}

//---------------------
// PUBLIC Clear vectors after filling trees
//---------------------
void Reconstruction::ClearTreeVectors() {

  // b-jet and b-quark vectors
  v_bPartons_Lorentz.clear();
  v_bPartons_Origin.clear();
  v_bPartons_pt.clear();
  v_bPartons_eta.clear();
  v_bPartons_phi.clear();
  v_bPartons_E.clear();
  v_bPartons_M.clear();

  // dR between b-quarks
  deltaR_bquarks.clear();

  v_SelectedJets_Lorentz.clear();
  v_SelectedJets_pt.clear();
  v_SelectedJets_eta.clear();
  v_SelectedJets_phi.clear();
  v_SelectedJets_E.clear();
  v_SelectedJets_M.clear();
  v_SelectedJets_bTag.clear();

  //results of jet-quark matching
  v_bIndexMatchedToJet.clear();
  v_JetMatchedOrigin.clear();

  return;
}


//---------------------
// PUBLIC Assign hard scatter particles to private vectors
//---------------------
void Reconstruction::SetHardScatterParticles(std::map<std::string,GenParticle*> in_GenParticles, std::map<std::string,TLorentzVector> in_GenParticlesLorentz, std::vector<std::pair<std::string,GenParticle*> > in_bQuarks, std::vector<std::pair<std::string,TLorentzVector> > in_bQuarksLorentz) {
  v_GenParticles = in_GenParticles;
  v_GenParticlesLorentz = in_GenParticlesLorentz;
  v_bQuarks = in_bQuarks;
  v_bQuarksLorentz = in_bQuarksLorentz;
  return;
}

//---------------------
// PUBLIC Assign leptons to private vectors
//---------------------
void Reconstruction::SetRecoLeptons(std::vector<std::pair<std::pair<std::string,int>, TLorentzVector> > in_RecoLeptonsLorentz) {
  v_RecoLeptonsLorentz = in_RecoLeptonsLorentz;
  return;
}

//---------------------
// PUBLIC Assign missing et to private MissingET*
//---------------------
void Reconstruction::SetMET(MissingET* in_MissingET) {
  m_MissingET = in_MissingET;
  return;
}

//---------------------
// PUBLIC Get reco jets and store those passing cuts
//---------------------
void Reconstruction::SetBJets(std::vector<Jet*> in_AllJets, std::vector<TLorentzVector> in_AllJetsLorentz, int num_btags, int num_maxbtags) {
  //  std::cout << "SetBJets()" << std::endl;
  n_btags = 0;

  // want event where there are at least 4 jets
  if ( in_AllJets.size() >= 4 && num_btags >= num_maxbtags) {
    v_AllJets = in_AllJets;
    v_AllJetsLorentz = in_AllJetsLorentz;
    std::sort(v_AllJets.begin(), v_AllJets.end(), SortByPt<Jet*>);
    std::vector<std::pair<int,Jet*> > v_btaggedJets;
    std::vector<std::pair<int,Jet*> > v_nonBTaggedJets;
    // get b-tagged jets
    std::vector<Jet*>::iterator jet_itr = v_AllJets.begin();
    std::vector<Jet*>::iterator jet_end = v_AllJets.end();
    for ( ; jet_itr != jet_end; ++jet_itr ) {
      int i_j = std::distance(v_AllJets.begin(), jet_itr);
      if ((*jet_itr)->BTag == 1) v_btaggedJets.push_back(std::make_pair(i_j, *jet_itr));
      else v_nonBTaggedJets.push_back(std::make_pair(i_j, *jet_itr));
    }

    std::vector<std::pair<int,Jet*> >::iterator bjet_itr = v_btaggedJets.begin();
    std::vector<std::pair<int,Jet*> >::iterator bjet_end = v_btaggedJets.end();
    for ( ; bjet_itr != bjet_end; ++bjet_itr ) {
      if (v_SelectedJets.size() < 4 ) {
        v_SelectedJets.push_back(*bjet_itr);
        v_SelectedJetsLorentz.push_back(std::make_pair((bjet_itr->first),(bjet_itr->second)->P4()));
        n_btags++;
      }
    }

    if (v_SelectedJets.size() < 4 ) {
      std::vector<std::pair<int,Jet*> >::iterator nonbjet_itr = v_nonBTaggedJets.begin();
      std::vector<std::pair<int,Jet*> >::iterator nonbjet_end = v_nonBTaggedJets.end();
      for ( ; nonbjet_itr != nonbjet_end; ++nonbjet_itr ) {
        if (v_SelectedJets.size() < 4) {
          v_SelectedJets.push_back(*nonbjet_itr);
          v_SelectedJetsLorentz.push_back(std::make_pair((nonbjet_itr->first),(nonbjet_itr->second)->P4()));
        }
      }
    }
    v_btaggedJets.clear();
    v_nonBTaggedJets.clear();
  }

  std::sort(v_SelectedJets.begin(),v_SelectedJets.end(), SortByPt2<int,Jet*>);
  std::sort(v_SelectedJetsLorentz.begin(),v_SelectedJetsLorentz.end(), SortLorentzByPt);

  return;
}


TTree* Reconstruction::GetTree() {
  return tree;
}

void Reconstruction::DeleteTrees() {
  delete tree;
  return;
}
