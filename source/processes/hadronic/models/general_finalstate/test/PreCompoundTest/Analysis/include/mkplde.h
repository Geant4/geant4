#ifndef mkplde_h
#define mkplde_h

#include "mkplsingledifhistogram.h"

using namespace std;

class mkplde : public mkplsingledifhistogram
{
private:
  inline mkplde();
public:
  inline mkplde(const int A, const int Z, const int pA, const int pZ,
		const int tA, const int tZ, const double E, const char * name);
  inline virtual ~mkplde();


  virtual void PrepareHistograms(mkplexpdatamanager * expdata);
  virtual void Normalization(const double xs, const double entries);
  virtual void Fill(const double E, const double ang, const string& pname);
  virtual vector<mkplcomparisonagreement> Plot(bool comp, const int n = -1);
  virtual void SaveEPS(const TString& bn, bool comp);

private:
  static const int mkpl_num_bins;

};

#include "mkplde.icc"

#endif
