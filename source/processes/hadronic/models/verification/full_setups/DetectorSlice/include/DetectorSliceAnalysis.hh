#ifndef DetectorSliceAnalysis_h 
#define DetectorSliceAnalysis_h 1

#include "globals.hh"
#include <vector>
#include <map>
#include <set>
#include "G4ThreeVector.hh"

#ifdef G4ANALYSIS_USE
namespace AIDA {
  class IAnalysisFactory;
  class IHistogramFactory;
  class ITree;
  class IHistogram1D;
  class IHistogram2D;
  class ITuple;
  class IPlotter;
}
#endif

class G4Step;
class G4Track;
class G4ParticleDefinition;


class DetectorSliceAnalysis {

public:

  static DetectorSliceAnalysis* getInstance();

  ~DetectorSliceAnalysis();
  
public:
  
  void init();
  void finish();
  void close();
  // The first and second methods are called by DetectorSliceRunAction, 
  // at the beginning and at the end of a Run, respectively.
  // The third method, close(), it is called at the end of
  // job, i.e. in the destructor of DetectorSliceRunAction.

  void fillNtuple( float incidentParticleId, 
		   float incidenteParticleEnergy, 
                   float totalEnergyDepositedInTracker,
                   float totalEnergyDepositedInActiveEmCal,
                   float totalEnergyDepositedInEmCal, 
                   float totalEnergyDepositedInActiveHadCal,
                   float totalEnergyDepositedInHadCal, 
                   float totalEnergyDepositedInMuonDet,
                   float exitingRadiusPrimaryMuon );
  // This method is supposed to be called at most once for each event.  

  void particleDepositingEnergy( const G4String caloType,
				 const G4double edep, const G4int particlePDG );
  // This method is called by DetectorSliceSensitiveCalorimeter
  // at each step in the active layer of the calorimeter.

  void infoTrack( const G4Track* aTrack );
  // This method is called by DetectorSliceTrackingAction
  // when each track is created. It is needed only for
  // the histograms of spectra of some particles types
  // when created in the two calorimeters.
 
  void endOfEvent( const G4double timeEventInSec );
  // Inform that the event if finished: this is useful to calculate
  // properly the statistical error of the quantities defined per
  // event.

private:
  
  DetectorSliceAnalysis();

  static DetectorSliceAnalysis* instance;
  
  G4int numberOfEvents;

  // To print, at the end of the Run, the average energy deposits, <E>,
  // and its error, in the whole calorimeter, in each readout layer, 
  // and in each radius bin, we need to collect the sum of the energy 
  // deposits and the sum of the square of energy deposits.
  // We also consider the maximum total energy deposit in an event,
  // to check for energy non conservation.
  G4int primaryParticleId;
  G4double beamEnergy;
  G4double sumEdepTracker, sumEdepTracker2;  // Deposited energy in the Tracker.
  G4double sumEdepEmAct, sumEdepEmAct2;      // Deposited energy in Sensitive EM Calo.
  G4double sumEdepEmTot, sumEdepEmTot2;      // Deposited energy in the all EM Calo.
  G4double sumEdepHadAct, sumEdepHadAct2;    // Deposited energy in Sensitive HAD Calo.
  G4double sumEdepHadTot, sumEdepHadTot2;    // Deposited energy in the all HAD Calo.
  G4double sumEdepMuon, sumEdepMuon2;        // Deposited energy in the Muon Detector.

  // Keep also separate information on the visible energy,
  // in EM and HAD calorimeters, for the following particles:
  G4double sumEdepEmAct_electron, sumEdepEmAct_electron2;    // e-   &  e+
  G4double sumEdepEmAct_muon, sumEdepEmAct_muon2;            // mu-  &  mu+
  G4double sumEdepEmAct_pion, sumEdepEmAct_pion2;            // pi-  &  pi+
  G4double sumEdepEmAct_kaon, sumEdepEmAct_kaon2;            // k-   &  k+
  G4double sumEdepEmAct_proton, sumEdepEmAct_proton2;        // p    &  pbar
  G4double sumEdepEmAct_nuclei, sumEdepEmAct_nuclei2;        // n    &  nuclei

  G4double sumEdepHadAct_electron, sumEdepHadAct_electron2;  // e-   &  e+
  G4double sumEdepHadAct_muon, sumEdepHadAct_muon2;          // mu-  &  mu+
  G4double sumEdepHadAct_pion, sumEdepHadAct_pion2;          // pi-  &  pi+
  G4double sumEdepHadAct_kaon, sumEdepHadAct_kaon2;          // k-   &  k+
  G4double sumEdepHadAct_proton, sumEdepHadAct_proton2;      // p    &  pbar
  G4double sumEdepHadAct_nuclei, sumEdepHadAct_nuclei2;      // n    &  nuclei

  // Due to the non-gaussian visible energy distribution, the width 
  // of it cannot be estimated as the sigma of the gaussian fit; on
  // the other hand, the rms is not appropriate because it weights 
  // too much the tails. A better estimate would be to consider the
  // "normalized deviation". To calculate it, we need to keep all
  // the values of visible energies in vectors.
  std::vector< G4double > vecEvisEm;
  std::vector< G4double > vecEvisHad;

  // To improve performance, save the PDG code of particles on 
  // constant integers. It is important that they are declared
  // as constant in order to use them in switch statements.
  const G4int electronId, positronId, gammaId;
  const G4int muonMinusId, muonPlusId, tauMinusId, tauPlusId;
  const G4int eNeutrinoId, antiENeutrinoId, muNeutrinoId, antiMuNeutrinoId,
    tauNeutrinoId, antiTauNeutrinoId;
  const G4int pionPlusId, pionMinusId, pionZeroId;
  const G4int kaonMinusId, kaonPlusId, kaonZeroId, antiKaonZeroId,
    kaonShortId, kaonLongId;
  const G4int protonId, antiProtonId, neutronId, antiNeutronId;

  // Monitor the CPU time event by event.
  std::multiset<G4double> eventTimeSet;

#ifdef G4ANALYSIS_USE
  AIDA::IAnalysisFactory* analysisFactory;
  AIDA::ITree* tree;
  AIDA::ITuple* tuple;
  AIDA::IHistogramFactory* histoFactory;

  // Kinetic spectra of some particles when created in EM and HAD calorimeters.
  AIDA::IHistogram1D* gammaSpectrumInEmCal;
  AIDA::IHistogram1D* neutronSpectrumInEmCal;
  AIDA::IHistogram1D* protonSpectrumInEmCal;
  AIDA::IHistogram1D* pionZeroSpectrumInEmCal;
  AIDA::IHistogram1D* pionPlusSpectrumInEmCal;
  AIDA::IHistogram1D* pionMinusSpectrumInEmCal;
                   
  AIDA::IHistogram1D* gammaSpectrumInHadCal;
  AIDA::IHistogram1D* neutronSpectrumInHadCal;
  AIDA::IHistogram1D* protonSpectrumInHadCal;
  AIDA::IHistogram1D* pionZeroSpectrumInHadCal;
  AIDA::IHistogram1D* pionPlusSpectrumInHadCal;
  AIDA::IHistogram1D* pionMinusSpectrumInHadCal;

#endif

};


#endif
