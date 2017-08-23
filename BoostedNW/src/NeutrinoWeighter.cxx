#include "include/NeutrinoWeighter.h"
#include <iostream>
#include "TObject.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

void NeutrinoWeighter::Reset(){

  m_top   = TLorentzVector();
  m_tbar  = TLorentzVector();
  m_ttbar = TLorentzVector();
  m_nu    = TLorentzVector();
  m_nubar = TLorentzVector();
  m_b     = TLorentzVector();
  m_bbar  = TLorentzVector();
  m_weight_max = -99.;

}

void NeutrinoWeighter::RecalculateEtas(double pos_lep_eta, double neg_lep_eta){

  m_nu_eta.clear();
  m_nu_sinh.clear();
  m_nu_cosh.clear();

  m_nubar_eta.clear();
  m_nubar_sinh.clear();
  m_nubar_cosh.clear();

  //UInt_t seed = 105200;
  //m_random.SetSeed(seed);

  for (int i = 0; i < 200; ++i){
    double mean_nu = pos_lep_eta*0.7 - 0.1;
    double eta_nu = m_random.Gaus(mean_nu, 1.14);

    m_nu_eta.push_back(eta_nu);
    m_nu_sinh.push_back(sinh(eta_nu));
    m_nu_cosh.push_back(cosh(eta_nu));

    double mean_nubar = neg_lep_eta*0.7 - 0.1;
    double eta_nubar = m_random.Gaus(mean_nubar, 1.14);

    m_nubar_eta.push_back(eta_nubar);
    m_nubar_sinh.push_back(sinh(eta_nubar));
    m_nubar_cosh.push_back(cosh(eta_nubar));

  }
}

NeutrinoWeighter::NeutrinoWeighter(){

  //std::cout << "Initialising Class" << std::endl;

  m_top   = TLorentzVector();
  m_tbar  = TLorentzVector();
  m_ttbar = TLorentzVector();
  m_b     = TLorentzVector();
  m_bbar  = TLorentzVector();
  m_nu    = TLorentzVector();
  m_nubar = TLorentzVector();
  m_random.SetSeed(105200);
  m_weight_max = -99.;
  m_do_both_pairings = false;

  ///-- Determine eta sampling --///
  for (double eta = -5.0; eta < 5.0001; eta += 0.2) {

    double sinh_eta = sinh(eta);
    double cosh_eta = cosh(eta);

    //std::cout << eta << " " << sinh_eta << " " << cosh_eta << std::endl;
    m_nu_eta.push_back(eta);
    m_nu_sinh.push_back(sinh_eta);
    m_nu_cosh.push_back(cosh_eta);

    m_nubar_eta.push_back(eta);
    m_nubar_sinh.push_back(sinh_eta);
    m_nubar_cosh.push_back(cosh_eta);
  }

  ///-- Do Top Mass Smearing --///
  //m_top_smears.push_back(169.0);
  //m_top_smears.push_back(169.5);
  //m_top_smears.push_back(170.0);
  //m_top_smears.push_back(170.5);
  //m_top_smears.push_back(170.75);
  //m_top_smears.push_back(171.0);
  //m_top_smears.push_back(171.25);
  //m_top_smears.push_back(171.5);
  //m_top_smears.push_back(171.75);
  m_top_smears.push_back(172.0);
  m_top_smears.push_back(172.25);
  m_top_smears.push_back(172.5);
  m_top_smears.push_back(172.75);
  m_top_smears.push_back(173.0);
  m_top_smears.push_back(173.25);
  m_top_smears.push_back(173.5);
  m_top_smears.push_back(173.75);
  m_top_smears.push_back(174.0);
  //m_top_smears.push_back(174.25);
  //m_top_smears.push_back(174.5);
  //m_top_smears.push_back(174.75);
  //m_top_smears.push_back(175.0);
  //m_top_smears.push_back(175.5);
  //m_top_smears.push_back(176.0);

  //m_W_smears.push_back(80.3);
  //m_W_smears.push_back(80.35);
  //m_W_smears.push_back(80.37);
  m_W_smears.push_back(80.4);
  //m_W_smears.push_back(80.45);
  //m_W_smears.push_back(80.5);


  //std::cout << "Random Seed = " << m_random.GetSeed() << std::endl;

}

double NeutrinoWeighter::Reconstruct(TLorentzVector lepton_pos, TLorentzVector lepton_neg, TLorentzVector jet_1, TLorentzVector jet_2, double met_ex, double met_ey, double met_phi){

  //m_mt2 = mt2(jet_1, jet_2, leptons[0], leptons[1], m_met_et, 30);
  double width_1, width_2;

  width_1 = 0.15;
  width_2 = 0.15;

  for(int mtop_counter = 0; mtop_counter < m_top_smears.size(); ++mtop_counter){
    for(int mtbar_counter = 0; mtbar_counter < m_top_smears.size(); ++mtbar_counter){
      for(int mWp_counter = 0; mWp_counter < m_W_smears.size(); ++mWp_counter){
        for(int mWn_counter = 0; mWn_counter < m_W_smears.size(); ++mWn_counter){

          double met_ex_original = met_ex;
          double met_ey_original = met_ey;

          calculateWeight(lepton_pos, lepton_neg, jet_1, jet_2, met_ex, met_ey, met_phi, m_top_smears[mtop_counter], m_top_smears[mtbar_counter], m_W_smears[mWp_counter], m_W_smears[mWn_counter]);
          if (m_do_both_pairings) {
           calculateWeight(lepton_pos, lepton_neg, jet_2, jet_1, met_ex, met_ey, met_phi, m_top_smears[mtop_counter], m_top_smears[mtbar_counter], m_W_smears[mWp_counter], m_W_smears[mWn_counter]);
          }

          ///-- Jet Smearing --///
          for (int smears = 0; smears < 50; ++smears){
            TLorentzVector jet_1_smeared, jet_2_smeared;

            jet_1_smeared = jet_1;
            jet_2_smeared = jet_2;

            //UInt_t seed = 105200;
            //m_random.SetSeed(seed);

            double scale_1 = m_random.Gaus(jet_1.Pt(), width_1*jet_1.Pt())/jet_1.Pt();
            double scale_2 = m_random.Gaus(jet_2.Pt(), width_2*jet_2.Pt())/jet_2.Pt();

            jet_1_smeared = jet_1_smeared*scale_1;
            jet_2_smeared = jet_2_smeared*scale_2;

            met_ex = met_ex_original + jet_1.Px() - jet_1_smeared.Px() + jet_2.Px() - jet_2_smeared.Px();
            met_ey = met_ey_original + jet_1.Py() - jet_1_smeared.Py() + jet_2.Py() - jet_2_smeared.Py();

            if (jet_1_smeared.M() > 0.0 && jet_2_smeared.M() > 0.0) {
              calculateWeight(lepton_pos, lepton_neg, jet_1_smeared, jet_2_smeared, met_ex, met_ey, met_phi, m_top_smears[mtop_counter], m_top_smears[mtbar_counter], m_W_smears[mWp_counter], m_W_smears[mWn_counter]);
              if (m_do_both_pairings) {
                calculateWeight(lepton_pos, lepton_neg, jet_2_smeared, jet_1_smeared, met_ex, met_ey, met_phi, m_top_smears[mtop_counter], m_top_smears[mtbar_counter], m_W_smears[mWp_counter], m_W_smears[mWp_counter]);
              }
            }

            met_ex = met_ex_original;
      	    met_ey = met_ey_original;
      	  } // END OF JET SMEARING
	      }// END OF Wn MASS SMEARING
      }// END OF Wp MASS SMEARING
    }// END OF TBAR MASS SMEARING
  }// END OF TOP MASS SMEARING

  return m_weight_max;
}

void NeutrinoWeighter::calculateWeight(TLorentzVector lepton_pos, TLorentzVector lepton_neg, TLorentzVector b1, TLorentzVector b2, double met_ex, double met_ey, double met_phi, double mtop, double mtbar, double mWp, double mWn) {

  double weight = -1.;
  TLorentzVector top, tbar, ttbar, Wp, Wn;

  //std::cout <<"Calc Weight " <<  mtop << " " << mtbar << " " << mWp << " " << mWn << std::endl;

  for (unsigned int nu_eta_index = 0; nu_eta_index < m_nu_eta.size(); ++nu_eta_index){

    //std::cout << "nu eta index " << nu_eta_index << std::endl;

    std::vector<TLorentzVector> solution_1 = solveForNeutrinoEta(&lepton_pos, &b1, nu_eta_index, 1, mtop, mWp);
    if(solution_1.size() == 0){
      continue; // if no solutions in sol1 then continue before calculating sol2;
    }

    for (unsigned int nubar_eta_index = 0; nubar_eta_index < m_nubar_eta.size(); ++nubar_eta_index){

      //std::cout << "nubar eta index " << nubar_eta_index << std::endl;

      std::vector<TLorentzVector> solution_2 = solveForNeutrinoEta(&lepton_neg, &b2, nubar_eta_index, -1, mtbar, mWn);
      if (solution_2.size() == 0){
	      continue; // sol2 has no solutions, continue;
      }


      //// SOLUTION 0, 0 ////
      weight = neutrino_weight(solution_1.at(0), solution_2.at(0), met_ex, met_ey, met_phi);
      //std::cout << "weight 0,0" << weight << std::endl;

      if(weight > m_weight_max  && weight > 0.000001){

        top   = (lepton_pos + b1 + solution_1.at(0));
      	tbar  = (lepton_neg + b2 + solution_2.at(0));

      	ttbar = top + tbar;

	///-- No point saving non-physical solutons, even if they have high weights --///
        if (ttbar.M() > 300. && top.E() > 0. && tbar.E() > 0.){
      	  m_weight_max = weight;
      	  m_top        = top;
      	  m_tbar       = tbar;
      	  m_ttbar      = ttbar;
      	  m_b          = b1;
      	  m_bbar       = b2;
      	  m_nu         = solution_1.at(0);
      	  m_nubar      = solution_2.at(0);
      	}
      }

      //// SOLUTION 1, 0 ////

      if (solution_1.size() > 1){

      	weight = neutrino_weight(solution_1.at(1), solution_2.at(0), met_ex, met_ey, met_phi);
      	//std::cout << "weight 1,0" << weight << std::endl;
      	if(weight > m_weight_max && weight > 0.000001){
      	  top   = (lepton_pos + b1 + solution_1.at(1));
      	  tbar  = (lepton_neg + b2 + solution_2.at(0));

      	  ttbar = top + tbar;

      	  ///-- No point saving non-physical solutons, even if they have high weights --///
       	  if (ttbar.M() > 300. && top.E() > 0 && tbar.E() > 0){
      	    m_weight_max = weight;
      	    m_top        = top;
      	    m_tbar       = tbar;
      	    m_ttbar      = ttbar;
      	    m_b          = b1;
      	    m_bbar       = b2;
      	    m_nu         = solution_1.at(1);
      	    m_nubar      = solution_2.at(0);
      	  }
      	}
      }

      //// SOLUTION 0, 1 ////
      if (solution_2.size() > 1){

      	weight = neutrino_weight(solution_1.at(0), solution_2.at(1), met_ex, met_ey, met_phi);
      	//std::cout << "weight 0,1" << weight << std::endl;
      	if(weight > m_weight_max && weight > 0.000001){
      	  top   = (lepton_pos + b1 + solution_1.at(0));
      	  tbar  = (lepton_neg + b2 + solution_2.at(1));

      	  ttbar = top + tbar;

      	  ///-- No point saving non-physical solutons, even if they have high weights --///
      	  if (ttbar.M() > 300. && top.E() > 0. && tbar.E() > 0.){
      	    m_weight_max = weight;
      	    m_top        = top;
      	    m_tbar       = tbar;
      	    m_ttbar      = ttbar;
      	    m_b          = b1;
      	    m_bbar       = b2;
      	    m_nu         = solution_1.at(0);
      	    m_nubar      = solution_2.at(1);
      	  }
      	}
      }

      //// SOLUTION 1, 1 ////
      if (solution_1.size() > 1 && solution_2.size() > 1){

      	weight = neutrino_weight(solution_1.at(1), solution_2.at(1), met_ex, met_ey, met_phi);
      	//std::cout << "weight 1,1" << weight << std::endl;
      	if(weight > m_weight_max && weight > 0.000001){
      	  top   = (lepton_pos + b1 + solution_1.at(1));
      	  tbar  = (lepton_neg + b2 + solution_2.at(1));

      	  ttbar = top + tbar;

      	  ///-- No point saving non-physical solutons, even if they have high weights --///
      	  if (ttbar.M() > 300. && top.E() > 0 && tbar.E() > 0){
      	    m_weight_max = weight;
      	    m_top        = top;
      	    m_tbar       = tbar;
      	    m_ttbar      = ttbar;
      	    m_b          = b1;
      	    m_bbar       = b2;
      	    m_nu         = solution_1.at(1);
      	    m_nubar      = solution_2.at(1);
      	  }
      	}
      }

    }//End of nubar eta index
  }// End of nu eta index

  return;
}

double NeutrinoWeighter::neutrino_weight(TLorentzVector neutrino1, TLorentzVector neutrino2, double met_ex, double met_ey, double met_phi){

  //if( abs(met_ex) < 0.01 || abs(met_ey) < 0.01)
  //  std::cout << "WARNING: one of the MET components is very very small! " << met_ex << met_ey << std::endl;

  double dx = met_ex - neutrino1.Px() - neutrino2.Px();
  double dy = met_ey - neutrino1.Py() - neutrino2.Py();
  double dphi = 1.;
  if(met_phi > -99.){
    TLorentzVector nunubar = neutrino1 + neutrino2;
    dphi = met_phi - nunubar.Phi();
  }


  //double m_sigma_met_ex = m_met_sumet*0.023 + 6.5;
  //double m_sigma_met_ey = m_met_sumet*0.023 + 6.5;
  double m_sigma_met_ex  = 0.2*met_ex;
  double m_sigma_met_ey  = 0.2*met_ey;
  double m_sigma_met_phi = 0.05;//0.15

  double numerator_x = -dx*dx;
  double numerator_y = -dy*dy;
  double numerator_phi = 1.;
  if(met_phi > -99.)
    numerator_phi = -dphi*dphi;

  double denominator_x = 2.*m_sigma_met_ex*m_sigma_met_ex;
  double denominator_y = 2.*m_sigma_met_ey*m_sigma_met_ey;
  double denominator_phi = 1.;
  if(met_phi > -99.)
    denominator_phi = 2.*m_sigma_met_phi*m_sigma_met_phi;

  double exp_x = exp(numerator_x/denominator_x);
  double exp_y = exp(numerator_y/denominator_y);
  double exp_phi = 1.;
  if(met_phi > 99.)
    exp_phi = exp(numerator_phi/denominator_phi);

  return exp_x*exp_y*exp_phi;

}

std::vector<TLorentzVector> NeutrinoWeighter::solveForNeutrinoEta(TLorentzVector* lepton, TLorentzVector* bJet, int index, int index_type, double mtop, double mW) {

  double nu_cosh = -99.;
  double nu_sinh = -99.;

  if (index_type > 0){
    nu_cosh = m_nu_cosh.at(index);
    nu_sinh = m_nu_sinh.at(index);
  } else {
    nu_cosh = m_nubar_cosh.at(index);
    nu_sinh = m_nubar_sinh.at(index);
  }

  //double Wmass2 = 80.4*80.4;
  double Wmass2 = mW*mW;
  double bmass = bJet->M();
  double Elprime = lepton->E() * nu_cosh - lepton->Pz() * nu_sinh;
  double Ebprime = bJet->E()   * nu_cosh - bJet->Pz()   * nu_sinh;

  double A = (lepton->Py() * Ebprime - bJet->Py() * Elprime) / (bJet->Px() * Elprime - lepton->Px() * Ebprime);
  double B = (Elprime * (mtop * mtop - Wmass2 - bmass * bmass - 2. * lepton->Dot(*bJet)) - Ebprime * Wmass2) / (2. * (lepton->Px() * Ebprime - bJet->Px() * Elprime));

  double par1 = (lepton->Px() * A + lepton->Py()) / Elprime;
  double C = A * A + 1. - par1 * par1;
  double par2 = (Wmass2 / 2. + lepton->Px() * B) / Elprime;
  double D = 2. * (A * B - par2 * par1);
  double F = B * B - par2 * par2;
  double det = D * D - 4. * C * F;


  std::vector<TLorentzVector> sol;

  ///-- 0 solutions case --///
  if (det < 0.0){
    return sol;                                                                                                                                                 }

  ///-- Only one real solution case --///
  if (det == 0.) {
    double py1 = -D / (2. * C);
    double px1 = A * py1 + B;
    double pT2_1 = px1 * px1 + py1 * py1;
    double pz1 = sqrt(pT2_1) * nu_sinh;

    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));

    if (!TMath::IsNaN(a1.E()) )
      sol.push_back(a1);
    return sol;
  }

  ///-- 2 solutions case --///
  if(det > 0){
    double tmp   = sqrt(det) / (2. * C);
    double py1   = -D / (2. * C) + tmp;
    double py2   = -D / (2. * C) - tmp;
    double px1   = A * py1 + B;
    double px2   = A * py2 + B;
    double pT2_1 = px1 * px1 + py1 * py1;
    double pT2_2 = px2 * px2 + py2 * py2;
    double pz1   = sqrt(pT2_1) * nu_sinh;
    double pz2   = sqrt(pT2_2) * nu_sinh;
    TLorentzVector a1(px1, py1, pz1, sqrt(pT2_1 + pz1 * pz1));
    TLorentzVector a2(px2, py2, pz2, sqrt(pT2_2 + pz2 * pz2));

    if (!TMath::IsNaN(a1.E()) && !TMath::IsNaN(a2.E())){
      sol.push_back(a1);
      sol.push_back(a2);
    }
    return sol;
  }

  ///-- Should never reach this point --///
  return sol;
}



/*double TopDileptonEventSaver::mt2(TLorentzVector jet_1, TLorentzVector jet_2, TLorentzVector lep_1, TLorentzVector lep_2, double met, int n_points){

  double mt2_result = 999999.;#

  TLorentzVector visible_1 = jet_1 + lep_1;
  TLorentzVector visible_2 = jet_2 + lep_2;

  double energy_1 = sqrt(pow(visible_1.M(), 2) + pow(visible_1.Pt(), 2));
  double energy_2 = sqrt(pow(visible_2.M(), 2) + pow(visible_2.Pt(), 2));

  for( int i = 0; i < n_points; ++i){

    double prob = (1./n_points)*i;
    double daughter_1 = prob*met;
    double daughter_2 = (1-prob)*met;

    double result = 99999999;
    double result_1 = sqrt( pow(visible_1.M(), 2) + 2*(energy_1*fabs(daughter_1) - (visible_1.Pt()*daughter_1)));
    double result_2 = sqrt( pow(visible_1.M(), 2) + 2*(energy_1*fabs(daughter_1) - (visible_1.Pt()*daughter_1)));

    if (result_1 > result_2) {
      result = result_1;
    } else {
      result = result_2;
    }

    if (result < mt2_result){
      mt2_result = result;
    }

  } /// End of probablility loop

  return mt2_result;

  }*/
