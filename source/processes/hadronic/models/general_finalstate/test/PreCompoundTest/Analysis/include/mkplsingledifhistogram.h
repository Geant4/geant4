#ifndef mkplsingledifhistogram_h
#define mkplsingledifhistogram_h 

#include "mkplvhistogram.h"

#include "TH1F.h"
#include "TMultiGraph.h"


using namespace std;

class mkplsingledifhistogram : public mkplvhistogram
{
protected:
  inline mkplsingledifhistogram();
public:
  inline mkplsingledifhistogram(const int A, const int Z, const int pA, const int pZ, 
				const int tA, const int tZ, const double E, 
				const char * name);
  inline virtual ~mkplsingledifhistogram();
  
  inline virtual void DeleteHistograms();
  inline virtual void InitializeHistograms();
  
  inline virtual bool ThereIsData() const;

  inline virtual void SaveHistograms();

protected:
  // Simulation histograms
  TH1F * mkpl_h;
  TH1F * mkpl_preeq_h;
  TH1F * mkpl_evap_h;
  TH1F * mkpl_fis_h;
  TH1F * mkpl_fermi_h;
  TH1F * mkpl_inc_h;
  // Experimental data
  TMultiGraph * mkpl_exp;
  // Weight
  double mkpl_weight;  
};

#include "mkplsingledifhistogram.icc"

#endif
