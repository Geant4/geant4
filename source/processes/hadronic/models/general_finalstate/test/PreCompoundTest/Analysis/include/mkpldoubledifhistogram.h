#ifndef mkpldoubledifhistogram_h
#define mkpldoubledifhistogram_h

#include <algorithm>

#include "mkplvhistogram.h"

#include "TH1F.h"
#include "TMultiGraph.h"


using namespace std;

class mkpldoubledifhistogram : public mkplvhistogram
{
protected:
  inline mkpldoubledifhistogram();
public:
  inline mkpldoubledifhistogram(const int A, const int Z, const int pA, const int pZ, 
				const int tA, const int tZ, const double E, 
				const char * name);
  inline virtual ~mkpldoubledifhistogram();

  inline virtual void DeleteHistograms();
  inline virtual void InitializeHistograms();

  inline int Size() const;

protected:

  // Simulation histograms
  vector<TH1F*> mkpl_h;
  vector<TH1F*> mkpl_preeq_h;
  vector<TH1F*> mkpl_evap_h;
  vector<TH1F*> mkpl_fis_h;
  vector<TH1F*> mkpl_fermi_h;
  vector<TH1F*> mkpl_inc_h;

  // Weight  
  vector<double> mkpl_weight;
  
  struct DeleteHistogramVector
  {
    template<typename T>
    void operator()(const T* ptr) const
    {
      delete ptr;
    }
  };
  
  
};

#include "mkpldoubledifhistogram.icc"

#endif
