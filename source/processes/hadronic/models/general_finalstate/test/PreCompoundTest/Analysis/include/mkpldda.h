#ifndef mkpldda_h
#define mkpldda_h

#include "mkpldoubledifhistogram.h"
#include "mkplexpdda.h"

using namespace std;

class mkpldda : public mkpldoubledifhistogram
{
private:
  inline mkpldda();
public:
  inline mkpldda(const int A, const int Z, const int pA, const int pZ,
		 const int tA, const int tZ, const double E, const char * name);
  inline virtual ~mkpldda();

  inline virtual void mkpldda::DeleteHistograms();
  inline virtual void mkpldda::InitializeHistograms();
  inline virtual bool mkpldda::ThereIsData() const;

  virtual void PrepareHistograms(mkplexpdatamanager * expdata);
  virtual void Normalization(const double xs, const double entries);
  virtual void Fill(const double E, const double ang, const string& pname);
  virtual vector<mkplcomparisonagreement> Plot(bool comp, const int n = -1);
  virtual void SaveHistograms();
  virtual void SaveEPS(const TString& bn, bool comp);

  inline pair<double,double> GetRange(const int n) const;

private:
  static const int mkpl_num_bins;
  
  // Experimental data
  vector<mkplexpdda> * mkpl_exp;

};

#include "mkpldda.icc"

#endif
