#ifndef NUCL_WATCHER_H
#define NUCL_WATCHER_H
#include "pair.h"
#include "vector"
#include <math.h>

class NuclWatcher {

public:

NuclWatcher(double z, vector<double> expa, vector<double> expcs, 
                     vector<double> experr, bool check, bool nucl) : nuclz(z),  
		     checkable(check), nucleable(nucl)  {
  exper_as = expa;
  exper_cs = expcs;
  exper_err = experr;

};

void watch(double a, double z) {

  const double small = 0.001;
  if(fabs(z - nuclz) < small) {
    bool here = false;
    for(int i = 0; i < simulated_as.size(); i++) {
      if(fabs(simulated_as[i] - a) < small) {
        simulated_cs[i] += 1.;
	here = true;
        break;
      };
    };
    if(!here) { simulated_as.push_back(a); simulated_cs.push_back(1.); };
  };
};

void setInuclCs(double csec, int nev) { 
  for(int i = 0; i < simulated_as.size(); i++) {
    double err = sqrt(simulated_cs[i])/simulated_cs[i];
    simulated_prob.push_back(simulated_cs[i]/nev);
    simulated_cs[i] *= csec/nev;
    simulated_errors.push_back(simulated_cs[i]*err);    
  };
};

double getChsq() const { return izotop_chsq; };

pair<double,double> getAverageRatio() const { 
 return pair<double,double>(average_ratio,aver_rat_err); };

pair<double,double> getExpCs() const {
  double cs = 0.;
  double err = 0.;
  for(int iz = 0; iz < exper_as.size(); iz++) {
    cs += exper_cs[iz];
    err += exper_err[iz];
  };
  return pair<double,double>(cs,err);
};

bool to_check() const { return checkable; };

bool look_forNuclei() const { return nucleable; };

pair<double,double> getInuclCs() const {
  double cs = 0.;
  double err = 0.;
  for(int iz = 0; iz < simulated_as.size(); iz++) {
    cs += simulated_cs[iz];
    err += simulated_errors[iz];
  };
  return pair<double,double>(cs,err);
};

void print() {
  const double small = 0.001;
  cout << endl << " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
  << endl;
  cout << " **** izotop Z **** " << nuclz << endl;
  izotop_chsq = 0.;
  average_ratio = 0.;
  aver_rat_err = 0.;
  double exp_cs = 0.;
  double exp_cs_err = 0.;
  double inucl_cs = 0.;
  double inucl_cs_err = 0.;
  vector<bool> not_used(simulated_cs.size(),true);
  int nmatched = exper_as.size();
  int nused = simulated_cs.size();
  double lhood = 0.;
  
  for(int iz = 0; iz < exper_as.size(); iz++) {
    double a = exper_as[iz];
    exp_cs += exper_cs[iz];
    exp_cs_err += exper_err[iz];
    bool found = false;
    for(int i = 0; i < simulated_as.size(); i++) {
      if(fabs(simulated_as[i] - a) < small) {
        double rat = simulated_cs[i]/exper_cs[iz];
	lhood += log10(rat)*log10(rat);
	double rat_err = sqrt(simulated_errors[i]*simulated_errors[i] +
	 exper_err[iz]*exper_err[iz]*rat*rat)/exper_cs[iz];
	average_ratio += rat;
	aver_rat_err += rat_err; 
        cout << " A " << a << " exp.cs " << exper_cs[iz] << " err " << 
	  exper_err[iz] << endl << 
	" sim. cs " << simulated_cs[i] << " err " << simulated_errors[i] << endl
	<< " ratio " << rat << " err " << rat_err << endl;
	cout << " simulated production rate " << simulated_prob[i] << endl;  	  
	not_used[i] = false;
	izotop_chsq += (rat - 1.)*(rat - 1.)/rat_err/rat_err; 
	found = true;
	nused--;
        break;
      };
    };
    if(!found) {
      cout << " not found exper.: A " << a << " exp.cs " << exper_cs[iz] 
           << " err " << exper_err[iz] << endl;
    }
     else {
      nmatched--;
    };
  };
  cout << " not found in simulations " << nmatched << endl;
  cout << " not found in exper: " << nused << endl;
    for(int i = 0; i < simulated_as.size(); i++) {
      inucl_cs += simulated_cs[i];
      inucl_cs_err += simulated_errors[i];
      if(not_used[i]) 
        cout << " extra simul.: A " << simulated_as[i] <<
	" sim. cs " << simulated_cs[i] << " err " << simulated_errors[i] 
	 << endl;
	cout << " simulated production rate " << simulated_prob[i] << endl;  	  
    };
  int matched = exper_as.size() - nmatched;
  if(matched > 0) {
    aver_lhood = lhood;
    aver_matched = matched;    
    lhood = pow(10.,sqrt(lhood/matched));
    cout << " matched " << matched << " CHSQ " << sqrt(izotop_chsq)/matched
                          << endl
      << " raw chsq " << izotop_chsq << endl
      << " average ratio " << average_ratio/matched 
      << " err " << aver_rat_err/matched << endl 
      << " lhood " << lhood << endl;
  }
   else {
    izotop_chsq = 0.;
    aver_lhood = 0.;    
  };    
  cout << " exper. cs " << exp_cs << " err " << exp_cs_err << endl
       << " inucl. cs " << inucl_cs << " err " << inucl_cs_err << endl;
  cout <<  " ++++++++++++++++++++++++++++++++++++++++++++++++++++++++ "
  << endl;
};

double getLhood() const { return aver_lhood; };

double getNmatched() const { return aver_matched; };

private: 

double nuclz;
double izotop_chsq;
double average_ratio;
double aver_rat_err;
double aver_lhood;
double aver_matched;

vector<double> exper_as;
vector<double> exper_cs;
vector<double> exper_err;

vector<double> simulated_as;
vector<double> simulated_cs;
vector<double> simulated_errors;
vector<double> simulated_prob;

bool checkable;
bool nucleable;

};        

#endif // NUCL_WATCHER_H 
