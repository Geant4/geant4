#ifndef mkpldd_h
#define mkpldd_h

#include "mkpldoubledifhistogram.h"

using namespace std;

class mkpldd : public mkpldoubledifhistogram
{
private:
  inline mkpldd();
public:
  inline mkpldd(const int A, const int Z, const int pA, const int pZ,
		const int tA, const int tZ, const double E, const char * name);
  inline virtual ~mkpldd();

  inline virtual void mkpldd::DeleteHistograms();
  inline virtual void mkpldd::InitializeHistograms();
  inline virtual bool mkpldd::ThereIsData() const;


  virtual void PrepareHistograms(mkplexpdatamanager * expdata);
  virtual void Normalization(const double xs, const double entries);
  virtual void Fill(const double E, const double ang, const string& pname);
  virtual vector<mkplcomparisonagreement> Plot(bool comp, const int n = -1);
  virtual void SaveHistograms();
  virtual void SaveEPS(const TString& bn, bool comp);

  inline double GetAngle(const int n) const;

private:
  void PlotAll();

  inline void ClearTemp();

private:
  static const int mkpl_num_bins;
  static const double mkpl_angular_window;

  // Experimental data
  vector<pair<double,TMultiGraph*> > * mkpl_exp;


  // Temporary data needed to plot all angles in the same picture.
  vector<pair<double,TMultiGraph*> > exp;
  vector<TH1F*> sim;
};

#include "mkpldd.icc"

#endif
