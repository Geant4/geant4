#ifndef mkplcomparisonhistograms_h
#define mkplcomparisonhistograms_h

#include <vector>
#include <utility>
#include <algorithm>

#include "TH1F.h"
#include "TString.h"
#include "TCanvas.h"

#include "mkplexpdatamanager.h"
#include "mkplcomparisonagreement.h"
#include "mkplexpdda.h"
#include "mkplda.h"
#include "mkplde.h"
#include "mkpldd.h"
#include "mkpldda.h"

using namespace std;

class TApplication;


class mkplcomparisonhistograms
{
public:
  inline mkplcomparisonhistograms(const int A, const int Z, const char * name);

  inline ~mkplcomparisonhistograms();
   
  inline void SetProjectileA(const int A);
  inline void SetProjectileZ(const int Z);
  inline void SetTargetA(const int A);
  inline void SetTargetZ(const int Z);
  inline void SetReactionE(const double E);
  
  inline void PrepareHistograms(mkplexpdatamanager * expdata, const bool components,
				const double xs, const double entries);

  inline void Fill(const int A, const int Z, const double E, const double ang, 
		   const double E_cm, const double ang_cm, const string& pname);

  inline TString GetName() const;
  
  inline bool ThereIsDE() const;
  inline bool ThereIsDEcm() const;
  inline bool ThereIsDA() const;
  inline bool ThereIsDAcm() const;
  inline bool ThereIsDD() const;
  inline bool ThereIsDDcm() const;
  inline bool ThereIsDDA() const;
  inline bool ThereIsDDAcm() const;
  inline int GetNumOfAngles() const;
  inline int GetNumOfAnglesCM() const;
  inline int GetNumOfRanges() const;
  inline int GetNumOfRangesCM() const;
  inline double GetAngle(const int n) const;
  inline double GetAngleCM(const int n) const;
  inline pair<double,double> GetRange(const int n) const;
  inline pair<double,double> GetRangeCM(const int n) const;

  vector<mkplcomparisonagreement> PlotDE(bool comp);
  vector<mkplcomparisonagreement> PlotDEcm(bool comp);
  vector<mkplcomparisonagreement> PlotDA(bool comp);
  vector<mkplcomparisonagreement> PlotDAcm(bool comp);

  vector<mkplcomparisonagreement> PlotDD(const int n, bool comp);
  vector<mkplcomparisonagreement> PlotDDcm(const int n, bool comp);
  vector<mkplcomparisonagreement> PlotDDA(const int n, bool comp);
  vector<mkplcomparisonagreement> PlotDDAcm(const int n, bool comp);

  void SaveHistograms(TFile * f);
  void SaveEPS(const TString& bn, bool comp);

private:
  
  inline mkplcomparisonhistograms();
  inline mkplcomparisonhistograms(const mkplcomparisonhistograms& right);
  inline const mkplcomparisonhistograms& operator=(const mkplcomparisonhistograms& right);

  inline TCanvas * RenewCanvas(TCanvas * canvas);


  void DeleteHistograms();
  void InitializeHistograms();

  void PrepareDEhistograms(mkplexpdatamanager * expdata, const bool components);
  void PrepareDAhistograms(mkplexpdatamanager * expdata, const bool components);
  void PrepareDDhistograms(mkplexpdatamanager * expdata, const bool components);
  void PrepareDDAhistograms(mkplexpdatamanager * expdata, const bool components);

  void DENormalization(const double xs, const double entries);
  void DANormalization(const double xs, const double entries);
  void DDNormalization(const double xs, const double entries);
  void DDANormalization(const double xs, const double entries);

  void FillDE(const double E, const double ang, const string& pname);
  void FillDEcm(const double E, const double ang, const string& pname);
  void FillDA(const double E, const double ang, const string& pname);
  void FillDAcm(const double E, const double ang, const string& pname);
  void FillDD(const double E, const double ang, const string& pname);
  void FillDDcm(const double E, const double ang, const string& pname);
  void FillDDA(const double E, const double ang, const string& pname);
  void FillDDAcm(const double E, const double ang, const string& pname);

  
  // ------------- DATA MEMBERS---------------

private:
  
  // Identify the kind of particles to compare to
  int mkpl_A;
  int mkpl_Z;
  TString mkpl_name;

  // This is the reaction data
  int mkpl_projectile_A;
  int mkpl_projectile_Z;
  int mkpl_target_A;
  int mkpl_target_Z;
  double mkpl_reaction_E;

  // DE histograms
  //--------------
  // simulation
  mkplde * mkpl_de;
  mkplde * mkpl_cm_de;


  // DA histograms
  //--------------
  // simulation
  mkplda * mkpl_da;
  mkplda * mkpl_cm_da;

  // DD histograms
  //--------------
  // simulation
  mkpldd * mkpl_dd;
  mkpldd * mkpl_cm_dd;


  // DDA histograms
  //---------------
  // simulation
  mkpldda * mkpl_dda;
  mkpldda * mkpl_cm_dda;

};

#include "mkplcomparisonhistograms.icc"

#endif
