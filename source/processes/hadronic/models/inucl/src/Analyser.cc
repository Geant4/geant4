#include "Analyser.h"
#include <math.h>

void Analyser::setInelCsec(double csec, bool withn) {

inel_csec = csec; // mb
withNuclei = withn;
cout << " total inelastic " << inel_csec << endl;

}

void Analyser::handleWatcherStatistics() {
const double small = 1.e-10;
cout << " ===================================================== " << endl;
cout << " ********** Izotop analysis ****************** " << endl;
double fgr = 0.;
double averat = 0.;
double ave_err = 0.;
double gl_chsq = 0.;
double tot_exper = 0.;
double tot_exper_err = 0.;
double tot_inucl = 0.;
double tot_inucl_err = 0.;
double checked = 0.;

for(int iw = 0; iw < ana_watchers.size(); iw++) {
  ana_watchers[iw].setInuclCs(inel_csec,eventNumber);
  ana_watchers[iw].print();
  if(ana_watchers[iw].to_check()) {
    pair<double,double> rat_err = ana_watchers[iw].getAverageRatio();
    averat += rat_err.first;
    ave_err += rat_err.second;
    gl_chsq += ana_watchers[iw].getChsq();   
    pair<double,double> cs_err = ana_watchers[iw].getExpCs();
    tot_exper += cs_err.first;
    tot_exper_err += cs_err.second;
    pair<double,double> inucl_cs_err = ana_watchers[iw].getInuclCs();
    tot_inucl += inucl_cs_err.first;
    tot_inucl_err += inucl_cs_err.second;
    double iz_checked = ana_watchers[iw].getNmatched();
    if(iz_checked > 0.) {
      fgr += ana_watchers[iw].getLhood();
      checked += iz_checked;    
    };
  };
};

if(checked > 0.) {
  gl_chsq = sqrt(gl_chsq)/checked;
  averat /= checked;
  ave_err /= checked;
  fgr = pow(10.,sqrt(fgr/checked)); 
};
cout << " ===================================================== " << endl;
cout << " total exper c.s. " << tot_exper << " err " << tot_exper_err <<
 " tot inucl c.s. " << tot_inucl << " err " << tot_inucl_err << endl;
cout << " checked total " << checked << " lhood " << fgr << endl
     << " average ratio " << averat << " err " << ave_err << endl
     << " global chsq " << gl_chsq << endl;
  
}

