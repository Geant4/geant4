#include "G4Analyser.hh"
#include <math.h>

void G4Analyser::setInelCsec(G4double csec, 
			     G4bool withn) {

  inel_csec = csec; // mb
  withNuclei = withn;

  if (verboseLevel > 1) {
    G4cout << " total inelastic " << inel_csec << G4endl;
  }
}

void G4Analyser::handleWatcherStatistics() {

  const G4double small = 1.0e-10;

  if (verboseLevel > 1) {
    G4cout << " ===================================================== " << G4endl;
    G4cout << " ********** Izotop analysis ****************** " << G4endl;
  };

  G4double fgr = 0.0;
  G4double averat = 0.0;
  G4double ave_err = 0.0;
  G4double gl_chsq = 0.0;
  G4double tot_exper = 0.0;
  G4double tot_exper_err = 0.0;
  G4double tot_inucl = 0.0;
  G4double tot_inucl_err = 0.0;
  G4double checked = 0.0;

  for(G4int iw = 0; iw < ana_watchers.size(); iw++) {
    ana_watchers[iw].setInuclCs(inel_csec, eventNumber);
    ana_watchers[iw].print();

    if(ana_watchers[iw].to_check()) {
      pair<G4double, G4double> rat_err = ana_watchers[iw].getAverageRatio();

      averat += rat_err.first;
      ave_err += rat_err.second;
      gl_chsq += ana_watchers[iw].getChsq();   

      pair<G4double, G4double> cs_err = ana_watchers[iw].getExpCs();

      tot_exper += cs_err.first;
      tot_exper_err += cs_err.second;

      pair<G4double, G4double> inucl_cs_err = ana_watchers[iw].getInuclCs();

      tot_inucl += inucl_cs_err.first;
      tot_inucl_err += inucl_cs_err.second;

      G4double iz_checked = ana_watchers[iw].getNmatched();

      if(iz_checked > 0.0) {
	fgr += ana_watchers[iw].getLhood();
	checked += iz_checked;    
      };
    };
  };

  if(checked > 0.0) {
    gl_chsq = sqrt(gl_chsq) / checked;
    averat /= checked;
    ave_err /= checked;
    fgr = pow(10.0, sqrt(fgr / checked)); 
  };

  if (verboseLevel > 1) {
    G4cout << " ===================================================== " << G4endl;
    G4cout << " total exper c.s. " << tot_exper << " err " << tot_exper_err <<
      " tot inucl c.s. " << tot_inucl << " err " << tot_inucl_err << G4endl;
    G4cout << " checked total " << checked << " lhood " << fgr << G4endl
	   << " average ratio " << averat << " err " << ave_err << G4endl
	   << " global chsq " << gl_chsq << G4endl;
  }
}

