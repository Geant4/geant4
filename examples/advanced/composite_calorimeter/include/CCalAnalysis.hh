///////////////////////////////////////////////////////////////////////////////
// File: CCalAnalysis.hh
// Description: CCalAnalysis is a singleton class and interfaces all user
//              analysis code
///////////////////////////////////////////////////////////////////////////////
#ifndef CCalAnalysis_h 
#define CCalAnalysis_h 1

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

namespace AIDA {
  class IAnalysisFactory;
  class IHistogramFactory;
  class ITree;
  class IHistogram1D;
  class IHistogram2D;
  class ITuple;
  class IPlotter;
}


class CCalAnalysis {
public:
  virtual ~CCalAnalysis();
  
public:
  void BeginOfRun(G4int n);
  void EndOfRun(G4int n);
  void EndOfEvent(G4int flag);

  void Init();
  void Finish();

  int maxbin() {return numberOfTimeSlices;}

  void InsertEnergyHcal(float*);
  void InsertEnergyEcal(float*);
  void InsertEnergy(float v);
  void InsertLateralProfile(float*);
  void InsertTime(float*); 

  void setNtuple(float* hcalE, float* ecalE, float elab, float x, float y, 
		 float z, float edep, float edec, float edhc);

  static CCalAnalysis* getInstance();

private:
  CCalAnalysis();
private:
  static CCalAnalysis* instance;

  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;
  AIDA::ITuple* tuple;

  enum {numberOfTimeSlices = 40}; 

  AIDA::IHistogram1D* energy;
  AIDA::IHistogram1D* hcalE[28];           // 28 hadronic modules
  AIDA::IHistogram1D* ecalE[49];           // 49 crystal towers
  AIDA::IHistogram1D* timeHist[40];        // 40 nanoseconds time window
  AIDA::IHistogram1D* lateralProfile[70];  // 70 centimeters lateral window
                                     // (indeed 64 should be enough)
};


#endif

