#ifndef MyAnalysis_h 
#define MyAnalysis_h 1

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


class MyAnalysis {

public:

  ~MyAnalysis();
  
public:
  
  void init();
  void finish();
  // These methods are supposed to be called, respectively,
  // at the beginning and at the end of a Run, 
  
  void fillNtuple( float incidentParticleId, 
		   float incidenteParticleEnergy, 
                   float totalEnergyDepositedInActiveLayers,
                   float totalEnergyDepositedInCalorimeter );
  // This method is supposed to be called at most once for each event.  

  static MyAnalysis* getInstance();
  
private:
  
  MyAnalysis();
  
  static MyAnalysis* instance;
  
  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;
  AIDA::ITuple* tuple;

};

#endif

