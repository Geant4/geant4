#ifndef mkplvhistogram_h
#define mkplvhistogram_h 1

#include <vector>
#include <string>

#include "TString.h"

class mkplexpdatamanager;
class TH1F;
class TMultiGraph;

#include "mkplcomparisonagreement.h"

using namespace std;

class mkplvhistogram
{
protected:
  inline mkplvhistogram();
public:
  inline mkplvhistogram(const int A, const int Z, const int pA, const int pZ, 
			const int tA, const int tZ, const double E, const char * name);
  inline virtual ~mkplvhistogram();

  virtual void DeleteHistograms() = 0;
  virtual void InitializeHistograms() = 0;
  virtual void PrepareHistograms(mkplexpdatamanager * expdata) = 0;
  virtual void Normalization(const double xs, const double entries) = 0;
  virtual void Fill(const double E, const double ang, const string& pname) = 0;
  virtual vector<mkplcomparisonagreement> Plot(const bool components, const int n = -1) = 0;
  virtual void SaveHistograms() = 0;
  virtual void SaveEPS(const TString& bn, bool comp) = 0;

  inline bool IsCM() const;
  inline void SetCM(const bool op = true);

  virtual bool ThereIsData() const = 0;

  vector<mkplcomparisonagreement> Agreement(TMultiGraph * exp, TH1F * sim);
  double GetMinimumNonNull(TH1F * h);

private:

  bool   mkpl_CM;

protected:

  int mkpl_A;
  int mkpl_Z;
  int mkpl_projectile_A;
  int mkpl_projectile_Z;
  int mkpl_target_A;
  int mkpl_target_Z;
  double mkpl_reaction_E;
  TString mkpl_name;

  // Model names
  static const char * mkpl_preeq_model_name;
  static const char * mkpl_evap_model_name;
  static const char * mkpl_fermi_model_name;
  static const char * mkpl_fis_model_name;
  static const char * mkpl_inc_model_name;

  // Comodity
  static const double degtorad;
  static const double radtodeg;

};

#include "mkplvhistogram.icc"

#endif
