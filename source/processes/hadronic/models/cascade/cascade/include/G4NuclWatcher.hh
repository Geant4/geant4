#ifndef G4NUCL_WATCHER_HH
#define G4NUCL_WATCHER_HH


#ifndef GLOB
#include "globals.hh"
#endif

#include "pair.h"
#include "vector"
#include <math.h>

class G4NuclWatcher {

public:

  G4NuclWatcher(G4double z, 
		vector<G4double> expa, 
		vector<G4double> expcs, 
		vector<G4double> experr, 
		G4bool check, 
		G4bool nucl)
    : nuclz(z),
      checkable(check),
      nucleable(nucl)  {

    exper_as = expa;
    exper_cs = expcs;
    exper_err = experr;
  };

  void watch(G4double a, 
	     G4double z) {

    const G4double small = 0.001;

    if(fabs(z - nuclz) < small) {
      G4bool here = false;

      for(G4int i = 0; i < simulated_as.size(); i++) {

	if(fabs(simulated_as[i] - a) < small) {
	  simulated_cs[i] += 1.0;
	  here = true;
	  break;
	};
      };
      if(!here) { 
	simulated_as.push_back(a); 
	simulated_cs.push_back(1.0); 
      };
    };
  };

  void setInuclCs(G4double csec, 
		  G4int nev) { 

    for(G4int i = 0; i < simulated_as.size(); i++) {
      double err = sqrt(simulated_cs[i]) / simulated_cs[i];

      simulated_prob.push_back(simulated_cs[i] / nev);
      simulated_cs[i] *= csec / nev;
      simulated_errors.push_back(simulated_cs[i] * err);    
    };
  };

  G4double getChsq() const { 

    return izotop_chsq; 
  };

  pair<G4double, G4double> getAverageRatio() const { 

    return pair<G4double, G4double>(average_ratio, aver_rat_err); 
  };

  pair<G4double, G4double> getExpCs() const {

    G4double cs = 0.0;
    G4double err = 0.0;

    for(G4int iz = 0; iz < exper_as.size(); iz++) {
      cs += exper_cs[iz];
      err += exper_err[iz];
    };

    return pair<G4double, G4double>(cs, err);
  };

  G4bool to_check() const { 

    return checkable; 
  };

  G4bool look_forNuclei() const { 

    return nucleable; 
  };

  pair<G4double, G4double> getInuclCs() const {

    G4double cs = 0.0;
    G4double err = 0.0;

    for(G4int iz = 0; iz < simulated_as.size(); iz++) {
      cs += simulated_cs[iz];
      err += simulated_errors[iz];
    };

    return pair<G4double, G4double>(cs, err);
  };

  void print() {

    const G4double small = 0.001;

    G4cout << G4endl << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
	   << G4endl;
    G4cout << " **** izotop Z **** " << nuclz << G4endl;

    izotop_chsq = 0.0;
    average_ratio = 0.0;
    aver_rat_err = 0.0;
    G4double exp_cs = 0.0;
    G4double exp_cs_err = 0.0;
    G4double inucl_cs = 0.0;
    G4double inucl_cs_err = 0.0;
    vector<G4bool> not_used(simulated_cs.size(), true);
    G4int nmatched = exper_as.size();
    G4int nused = simulated_cs.size();
    G4double lhood = 0.0;
  
    for(G4int iz = 0; iz < exper_as.size(); iz++) {
      G4double a = exper_as[iz];

      exp_cs += exper_cs[iz];
      exp_cs_err += exper_err[iz];

      G4bool found = false;

      for(G4int i = 0; i < simulated_as.size(); i++) {

	if(fabs(simulated_as[i] - a) < small) {
	  G4double rat = simulated_cs[i] / exper_cs[iz];

	  lhood += log10(rat) * log10(rat);

	  G4double rat_err = sqrt(simulated_errors[i] * simulated_errors[i] +
				  exper_err[iz] * exper_err[iz] * rat * rat) / exper_cs[iz];
	  average_ratio += rat;
	  aver_rat_err += rat_err; 

	  G4cout << " A " << a << " exp.cs " << exper_cs[iz] << " err " << 
	    exper_err[iz] << G4endl << 
	    " sim. cs " << simulated_cs[i] << " err " << simulated_errors[i] << G4endl
		 << " ratio " << rat << " err " << rat_err << G4endl;
	  G4cout << " simulated production rate " << simulated_prob[i] << G4endl;  	  

	  not_used[i] = false;
	  izotop_chsq += (rat - 1.0) * (rat - 1.0) / rat_err / rat_err; 
	  found = true;
	  nused--;
	  break;
	};
      };
      if(!found) {
	G4cout << " not found exper.: A " << a << " exp.cs " << exper_cs[iz] 
	       << " err " << exper_err[iz] << G4endl;
      }
      else {
	nmatched--;
      };
    };

    G4cout << " not found in simulations " << nmatched << G4endl;
    G4cout << " not found in exper: " << nused << G4endl;

    for(G4int i = 0; i < simulated_as.size(); i++) {
      inucl_cs += simulated_cs[i];
      inucl_cs_err += simulated_errors[i];

      if(not_used[i]) 
        G4cout << " extra simul.: A " << simulated_as[i] <<
	  " sim. cs " << simulated_cs[i] << " err " << simulated_errors[i] << G4endl;
      G4cout << " simulated production rate " << simulated_prob[i] << G4endl;  	  
    };
    G4int matched = exper_as.size() - nmatched;

    if(matched > 0) {
      aver_lhood = lhood;
      aver_matched = matched;    
      lhood = pow(10.0, sqrt(lhood / matched));

      G4cout << " matched " << matched << " CHSQ " << sqrt(izotop_chsq) / matched
	     << G4endl
	     << " raw chsq " << izotop_chsq << G4endl
	     << " average ratio " << average_ratio / matched 
	     << " err " << aver_rat_err / matched << G4endl 
	     << " lhood " << lhood << G4endl;
    }
    else {
      izotop_chsq = 0.0;
      aver_lhood = 0.0;    
    }; 
   
    G4cout << " exper. cs " << exp_cs << " err " << exp_cs_err << G4endl
	   << " inucl. cs " << inucl_cs << " err " << inucl_cs_err << G4endl;
    G4cout <<  " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ " << G4endl;
  };

  G4double getLhood() const { 

    return aver_lhood; 
  };

  G4double getNmatched() const { 

    return aver_matched; 
  };

private: 

  G4double nuclz;

  G4double izotop_chsq;

  G4double average_ratio;

  G4double aver_rat_err;

  G4double aver_lhood;

  G4double aver_matched;

  vector<G4double> exper_as;

  vector<G4double> exper_cs;

  vector<G4double> exper_err;

  vector<G4double> simulated_as;

  vector<G4double> simulated_cs;

  vector<G4double> simulated_errors;

  vector<G4double> simulated_prob;

  G4bool checkable;

  G4bool nucleable;

};        

#endif // G4NUCL_WATCHER_HH 

