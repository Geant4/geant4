//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// 
// ------------------------------------------------------------
//        GEANT 4  include file implementation
// ------------------------------------------------------------
//
// Class description:
//
// G4CoupledTransportation is an optional process to transport  
// a particle, in case of coupled navigation in parallel geometries
//  i.e. the geometrical propagation will be done
//   encountering the geometrical volumes of the detectors and
//   those of parallel geometries (eg for biasing, scoring, fast simulation)
// It is tasked with updating the "safety" to reflect the geometrical
//   distance to the nearest volume, and the time of flight of the particle.

// =======================================================================
// Created:  17 May 2006, J. Apostolakis
// =======================================================================
#ifndef G4CoupledTransportation_hh
#define G4CoupledTransportation_hh 1

#include "G4VProcess.hh"

#include "G4FieldManager.hh"

#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4PathFinder.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleChangeForTransport.hh"

class G4SafetyHelper; 
class G4TransportationLogger;

class G4CoupledTransportation : public G4VProcess 
{
  // Concrete class that does the geometrical transport 

  public:  // with description

     G4CoupledTransportation( G4int verbosityLevel= 0); 
     ~G4CoupledTransportation(); 

     G4double      AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                                   G4double  previousStepSize,
                                   G4double  currentMinimumStep, 
                                   G4double& currentSafety,
                                   G4GPILSelection* selection
                            );

     G4VParticleChange* AlongStepDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            );

     G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step&  stepData
                            );
       // Responsible for the relocation

     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4double   previousStepSize,
                             G4ForceCondition* pForceCond
                            );
       // Forces the PostStepDoIt action to be called, 
       // but does not limit the step

     G4PropagatorInField* GetPropagatorInField();
     void SetPropagatorInField( G4PropagatorInField* pFieldPropagator);
       // Access/set the assistant class that Propagate in a Field

     inline G4double GetThresholdWarningEnergy() const; 
     inline G4double GetThresholdImportantEnergy() const; 
     inline G4int GetThresholdTrials() const; 

     inline void SetThresholdWarningEnergy( G4double newEnWarn ); 
     inline void SetThresholdImportantEnergy( G4double newEnImp ); 
     inline void SetThresholdTrials(G4int newMaxTrials ); 

     void SetHighLooperThresholds(); // Shortcut method - old values (meant for HEP)
     void SetLowLooperThresholds(); // Set low thresholds - for low-E applications   
     void PushThresholdsToLogger(); // Inform logger of current thresholds
     void ReportLooperThresholds(); // Print values of looper thresholds
   
     // Get/Set parameters for killing loopers: 
     //   Above 'important' energy a 'looping' particle in field will 
     //   *NOT* be abandoned, except after fThresholdTrials attempts.
     // Below Warning energy, no verbosity for looping particles is issued

     inline G4double GetMaxEnergyKilled() const; 
     inline G4double GetSumEnergyKilled() const;
     inline void ResetKilledStatistics( G4int report = 1);      
     // Statistics for tracks killed (currently due to looping in field)

     static G4bool EnableMagneticMoment(G4bool useMoment=true);    
     // Whether to deflect particles with force due to magnetic moment

     static G4bool EnableGravity(G4bool useGravity);
     // Turn on the capability to deflect particles with a gravity field

     static void   SetSilenceLooperWarnings( G4bool val);
     // Do not warn (or throw exception) about 'looping' particles
     static G4bool GetSilenceLooperWarnings();
   
     static void  SetSignifyStepsInAnyVolume( G4bool anyVol )
       { fSignifyStepInAnyVolume = anyVol; } 
     static G4bool GetSignifyStepsInAnyVolume()
       { return fSignifyStepInAnyVolume; }
     // Flag in step corresponds to first/last step in a volume 'any'
     // geometry (if this is true) or refers to first/last step in mass
     // geometry only (if false)

     // The following methods give access to first/last step in particular
     // geometry *independent* of the choice of the 'Signify' flag
     //
     G4bool IsFirstStepInAnyVolume() const { return fFirstStepInAnyVolume; }
     G4bool IsLastStepInAnyVolume() const { return fAnyGeometryLimitedStep; } 
     G4bool IsFirstStepInMassVolume() const { return fFirstStepInMassVolume; }
     G4bool IsLastStepInMassVolume() const { return fMassGeometryLimitedStep; } 

  public:  // without description

     void StartTracking(G4Track* aTrack); 
     void EndTracking();

     static G4bool EnableUseMagneticMoment(G4bool useMoment=true)
      { return EnableMagneticMoment(useMoment); }
     // Old name ... obsolete
   
     G4double AtRestGetPhysicalInteractionLength( const G4Track& ,
                                                  G4ForceCondition* )
       { return -1.0; }   // No operation in AtRestGPIL

     G4VParticleChange* AtRestDoIt( const G4Track&, const G4Step&)
       { return 0; }      // No operation in AtRestDoIt

     void PrintStatistics( std::ostream& outStr) const;
   
  protected:

     G4bool DoesAnyFieldExist();
       // Check whether any field exists in the geometry
       //  - replaces method that checked only whether a field for the world volume

     void ReportInexactEnergy(G4double startEnergy, G4double endEnergy);
       // Issue warning

     void ReportMove( G4ThreeVector OldVector, G4ThreeVector NewVector,
                      const G4String& Quantity );
   
     void ReportMissingLogger(const char * methodName);
   
  private:

     G4Navigator*         fMassNavigator;
       // The navigator for the 'mass' geometry
       // (the real one, that physics occurs in)
     G4PathFinder*        fPathFinder;
     G4int fNavigatorId;
       // The PathFinder used to transport the particle

     G4PropagatorInField* fFieldPropagator;
       // Still required in order to find/set the fieldmanager

     G4bool               fAnyFieldExists; 
     // G4bool fStartedNewTrack;   // True for first step or restarted tracking
                                   // until first step's AlongStepGPIL

     G4ThreeVector        fTransportEndPosition;
     G4ThreeVector        fTransportEndMomentumDir;
     G4double             fTransportEndKineticEnergy;
     G4ThreeVector        fTransportEndSpin;
     G4bool               fMomentumChanged;
       // The particle's state after this Step, Store for DoIt

     G4bool               fEndGlobalTimeComputed; 
     G4double             fCandidateEndGlobalTime;

     G4bool               fParticleIsLooping;
     G4bool               fNewTrack;
     G4ThreeVector        fPreviousSftOrigin; 
     G4double             fPreviousMassSafety;
     G4double             fPreviousFullSafety;

     G4TouchableHandle    fCurrentTouchableHandle;
     
     // G4bool         fFieldExists;
       // Whether a magnetic field exists ...
       // A data member for this is problematic: it is useful only if it
       // can be initialised and updated -- and a scheme is not yet possible.

     G4bool fMassGeometryLimitedStep;
       // Flag to determine whether a 'mass' boundary was reached.
     G4bool fAnyGeometryLimitedStep; 
       // Did any geometry limit the step ?

     G4ParticleChangeForTransport fParticleChange;
       // New ParticleChange

     G4double fEndpointDistance;


     // Thresholds for looping particles: 
     // 
     G4double fThreshold_Warning_Energy = 1.0 * CLHEP::keV;   //  Warn above this energy
     G4double fThreshold_Important_Energy = 1.0 * CLHEP::MeV; //  Give a few trial above this E
     G4int    fThresholdTrials = 10;       //  Number of trials an important looper survives   
       // Above 'important' energy a 'looping' particle in field will 
       // *NOT* be abandoned, except after fThresholdTrials attempts.

     // Counter for steps in which particle reports 'looping',
     // if it is above 'Important' Energy
     //
     G4int fNoLooperTrials=0;

     // Statistics for tracks abandoned
     //
     G4double fSumEnergyKilled= 0.0;
     G4double fSumEnerSqKilled= 0.0;
     G4double fMaxEnergyKilled= -1.0;
     G4int    fMaxEnergyKilledPDG= 0;
     unsigned long fNumLoopersKilled= 0;
     G4double fSumEnergyKilled_NonElectron= 0.0;
     G4double fSumEnerSqKilled_NonElectron= 0.0;
     G4double fMaxEnergyKilled_NonElectron= -1.0;
     G4int    fMaxEnergyKilled_NonElecPDG= 0;
     unsigned long fNumLoopersKilled_NonElectron= 0;
     G4double fSumEnergySaved=  0.0;
     G4double fMaxEnergySaved= -1.0;
     G4double fSumEnergyUnstableSaved = 0.0;
   
     G4SafetyHelper* fpSafetyHelper;  // To pass it the safety value obtained
     G4TransportationLogger* fpLogger;  // Reports issues / raises warnings
   
  private:

     friend class G4Transportation;
     static G4bool fUseMagneticMoment;
     static G4bool fUseGravity;
     static G4bool fSilenceLooperWarnings;  // Flag to *Supress* all 'looper' warnings   

  private:

     G4bool fFirstStepInMassVolume;
     G4bool fFirstStepInAnyVolume;
     // G4bool fLastStepInMassVolume; => use fMassGeometryLimitedStep 
     // G4bool fLastStepInAnyVolume; => use fAnyGeometryLimitedStep

     static G4bool fSignifyStepInAnyVolume;
       // True: First/Last step in any one of the geometries
       // False: First/Last step in volume of 'mass' geometry
};

#include "G4CoupledTransportation.icc"

#endif  
