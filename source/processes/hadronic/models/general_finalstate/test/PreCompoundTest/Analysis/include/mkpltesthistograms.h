#ifndef mkpltesthistograms_h
#define mkpltesthistograms_h

#include <vector>
#include <set>
#include <string>

#include "TVector3.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;

class mkpltesthistograms
{
public:
  inline mkpltesthistograms();
  inline ~mkpltesthistograms();

  void PrepareHistograms(const TVector3& IncidentDirection,
			 const set<string>& partt,
			 const set<string>& proct);

  void FillConservationHistograms(const int deltaA, 
				  const int deltaZ, 
				  const double deltaE, 
				  const double deltaP, 
				  const TVector3 & MomentumTest);

  void FillHistograms(const double Theta, 
		      const TVector3 & fragP, 
		      const string & procname,
		      const string & fname);

  void Draw(TApplication * app);

  inline void PlotBaryonNumber();
  inline void PlotCharge();
  inline void PlotEnergy();
  inline void PlotMomentum();
  inline void PlotPx();
  inline void PlotPy();
  inline void PlotPz();
  inline void PlotTheta();
  inline void PlotThetaPreeq();
  inline void PlotThetaEvap();
  inline void PlotPhi();
  inline void PlotPhiPreeq();
  inline void PlotPhiEvap();
  inline void PlotPhiNucleons();
  inline void PlotTypeOfFragments();
    
private:
  inline mkpltesthistograms(const mkpltesthistograms& right);
  inline const mkpltesthistograms& operator=(const mkpltesthistograms& right);

  void DeleteHistograms();
  void TestPhi(const TVector3& IncidentDirection);

  inline TCanvas * RenewCanvas(TCanvas * canvas);

private:
  bool mkpl_testPhi;
  char mkpl_direction;

  vector<string> mkpl_particletypes;
  vector<string> mkpl_processtypes;

  // Basic Test Histograms
  // ---------------------
  // Histogram for Baryonic number conservation
  TH1F * mkpl_baryonconservation;
  // Histogram for Charge conservation
  TH1F * mkpl_chargeconservation;
  // Histogram for Energy Conservation
  TH1F * mkpl_energyconservation;
  // Histograms for Momentum Conservation
  TH1F * mkpl_momentumconservation;
  TH1F * mkpl_pxconservation;
  TH1F * mkpl_pyconservation;
  TH1F * mkpl_pzconservation;
  // Histograms for Angular tests
  TH1F * mkpl_theta;
  TH1F * mkpl_thetaprecom;
  TH1F * mkpl_thetaevap;
  TH1F * mkpl_phi;
  TH1F * mkpl_phiprecom;
  TH1F * mkpl_phievap;
  TH1F * mkpl_phinucleon;
  
  // Histogram for kind of particles
  TH1F * mkpl_typefragments;

  // Constants for histograms
  static const int mkpl_num_bins_test;
};

#include "mkpltesthistograms.icc"

#endif
