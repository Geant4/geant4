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
	     const G4int numberOfReadoutLayersIn,
	     const G4int numberOfRadiusBinsIn,
	     const G4double radiusBinIn );
  void init();
  void finish();
  void close();
  // The first method is called by StatAccepTestDetectorConstruction
  // only when the geometry is constructed the first time or changed.  
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
  AIDA::IHistogramFactory* histoFactory;

  G4int numberOfEvents;

  G4int numberOfReplicas;
  G4int numberOfReadoutLayers;
  G4int numberOfActiveLayersPerReadoutLayer;
  G4int numberOfRadiusBins;  // Number of bins in the transverse profile
  G4double radiusBin;        // Size of the first bin in the transverse profile.
                             // The other bins have a bin size that is 
                             // as follows: second bin has a size that is twice
                             // the first bin size; the third one has a size
                             // that is three times the first bin size, and so on.

  std::vector< G4double > longitudinalProfile;
  std::vector< G4double > transverseProfile; 

  // To print, at the end of the Run, the average energy deposits, <E>,
  // and its error, in the whole calorimeter, in each readout layer, 
  // and in each radius bin, we need to collect the sum of the energy 
  // deposits and the sum of the square of energy deposits.
  // We also consider the maximum total energy deposit in an event,
  // to check for energy non conservation.
  G4double beamEnergy;
  G4double sumEdepAct, sumEdepAct2;
  G4double sumEdepTot, sumEdepTot2, maxEdepTot;
  std::vector< G4double > sumL;
  std::vector< G4double > sumL2;
  std::vector< G4double > sumR;
  std::vector< G4double > sumR2;

  // Due to the non-gaussian visible energy distribution, the width 
  // of it cannot be estimated as the sigma of the gaussian fit; on
  // the other hand, the rms is not appropriate because it weights 
  // too much the tails. A better estimate would be to consider the
  // "normalized deviation". To calculate it, we need to keep all
  // the values of EdepAct in a vector.
  std::vector< G4double > vecEvis;

  // Summary histograms for the longitudinal and transverse shower profiles. 
  AIDA::IHistogram1D* longitudinalProfileHisto;
  AIDA::IHistogram1D* transverseProfileHisto;

};

#endif

