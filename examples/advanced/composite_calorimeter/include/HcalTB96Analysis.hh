#ifndef HcalTB96Analysis_h 
#define HcalTB96Analysis_h 1

#include "globals.hh"
#include "g4std/vector"
#include "G4ThreeVector.hh"

class IAnalysisFactory;
class IHistogramFactory;
class ITree;
class IHistogram1D;
class IHistogram2D;
class ITuple;
class IPlotter;

class HcalTB96Analysis {
public:
  virtual ~HcalTB96Analysis();
  
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

  static HcalTB96Analysis* getInstance();

private:
  HcalTB96Analysis();
private:
  static HcalTB96Analysis* instance;

  IAnalysisFactory* analysisFactory;
  ITree* tree;
  ITuple* tuple;

  enum {numberOfTimeSlices = 40};

  IHistogram1D* energy;
  IHistogram1D* profile;
  IHistogram1D* hcalE[28];
  IHistogram1D* ecalE[49];
  IHistogram1D* timeHist[40];
  IHistogram1D* lateralProfile[28];

};


#endif

