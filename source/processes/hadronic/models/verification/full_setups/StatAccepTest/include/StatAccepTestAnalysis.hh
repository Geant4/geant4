#ifndef StatAccepTestAnalysis_h 
#define StatAccepTestAnalysis_h 1

#include "globals.hh"
#include <vector>
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


class StatAccepTestAnalysis {

public:

  ~StatAccepTestAnalysis();
  
public:
  
  void init( const G4int numberOfReplicasIn, 
	     const G4int numberOfRadiusBinsIn,
	     const G4double radiusBinIn );
  void init();
  void finish();
  void close();
  // The first method is called by StatAccepTestDetectorConstruction
  // only when the geometry is changed.  
  // The second and third methods are called by StatAccepTestRunAction, 
  // at the beginning and at the end of a Run, respectively.
  // The fourth method, close(), it is called at the end of
  // job, i.e. in the destructor of StatAccepTestRunAction.

  void fillNtuple( float incidentParticleId, 
		   float incidenteParticleEnergy, 
                   float totalEnergyDepositedInActiveLayers,
                   float totalEnergyDepositedInCalorimeter );
  // This method is supposed to be called at most once for each event.  

  void fillShowerProfile( G4int replica, G4double radius, G4double edep );
  // This method is called by StatAccepTestSensitiveCalorimeter at each step
  // in the active layer of the calorimeter.

  static StatAccepTestAnalysis* getInstance();

private:
  
  StatAccepTestAnalysis();
  
  static StatAccepTestAnalysis* instance;
  
  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;
  AIDA::ITuple* tuple;

  G4int numberOfReplicas;
  G4int numberOfRadiusBins;  // Number of bins in the transverse profile
  G4int numberOfEvents;
  G4double radiusBin;        // Size of the first bin in the transverse profile.
                             // The other bins have a bin size that is 
                             // as follows: second bin has a size that is twice
                             // the first bin size; the third one has a size
                             // that is three times the first bin size, and so on.

  std::vector< G4double > longitudinalProfile;
  std::vector< G4double > transverseProfile; 

  // To print, at the end of the Run, the average energy deposits, <E>,
  // and its error, in all active layers, in all calorimeter, in
  // each active layer, and in each radius bin, we need to collect 
  // the sum of the energy deposits and the sum of the square of 
  // energy deposits.
  G4double sumEdepAct, sumEdepAct2;
  G4double sumEdepTot, sumEdepTot2;
  std::vector< G4double > sumL;
  std::vector< G4double > sumL2;
  std::vector< G4double > sumR;
  std::vector< G4double > sumR2;

};

#endif

