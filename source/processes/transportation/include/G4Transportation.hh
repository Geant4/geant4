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
// G4Transportation is a process responsible for the transportation of 
// a particle, i.e. the geometrical propagation encountering the 
// geometrical sub-volumes of the detectors.
// It is also tasked with part of updating the "safety".

// =======================================================================
// Created:  19 March 1997, J. Apostolakis
// =======================================================================
#ifndef G4Transportation_hh
#define G4Transportation_hh 1

#include "G4VProcess.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleChangeForTransport.hh"

class G4Navigator;
class G4PropagatorInField;
class G4SafetyHelper; 
class G4TransportationLogger;

class G4Transportation : public G4VProcess 
{
  // Concrete class that does the geometrical transport 

  public:  // with description

     G4Transportation( G4int verbosityLevel= 1, const G4String& aName = "Transportation");
     ~G4Transportation(); 

     G4double      AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                                   G4double  previousStepSize,
                                   G4double  currentMinimumStep, 
                                   G4double& currentSafety,
                                   G4GPILSelection* selection
                            ); // override;

     G4VParticleChange* AlongStepDoIt(
                             const G4Track& track,
                             const G4Step& stepData
                            ); // override; 

     G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step&  stepData
                            ); // override;    
       // Responsible for the relocation

     G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& ,
                             G4double   previousStepSize,
                             G4ForceCondition* pForceCond
                            ); // override;                            
       // Forces the PostStepDoIt action to be called, 
       // but does not limit the step

     inline G4bool FieldExertedForce() { return fFieldExertedForce; }
   
     G4PropagatorInField* GetPropagatorInField();
     void SetPropagatorInField( G4PropagatorInField* pFieldPropagator);
       // Access/set the assistant class that Propagate in a Field

     inline G4double GetThresholdWarningEnergy() const; 
     inline G4double GetThresholdImportantEnergy() const; 
     inline G4int GetThresholdTrials() const; 

     inline void SetThresholdWarningEnergy( G4double newEnWarn ); 
     inline void SetThresholdImportantEnergy( G4double newEnImp ); 
     inline void SetThresholdTrials(G4int newMaxTrials ); 
     // Get/Set parameters for killing loopers: 
     //   Above 'important' energy a 'looping' particle in field will 
     //   *NOT* be abandoned, except after fThresholdTrials attempts.
     // Below Warning energy, no verbosity for looping particles is issued

     void SetHighLooperThresholds(); // Shortcut method - old values (meant for HEP)   
     void SetLowLooperThresholds(); // Set low thresholds - for low-E applications
     void PushThresholdsToLogger(); // Inform logger of current thresholds
     void ReportLooperThresholds(); // Print values of looper thresholds

     inline G4double GetMaxEnergyKilled() const; 
     inline G4double GetSumEnergyKilled() const;
     inline void ResetKilledStatistics( G4int report = 1);      
     // Statistics for tracks killed (currently due to looping in field)

     inline void EnableShortStepOptimisation(G4bool optimise=true); 
     // Whether short steps < safety will avoid to call Navigator (if field=0)

     static G4bool EnableMagneticMoment(G4bool useMoment=true); 
     // Whether to enable particles to be deflected with force due to magnetic moment

     static G4bool EnableGravity(G4bool useGravity=true); 
     // Whether to enable particles to be deflected with force due to gravity

     static void   SetSilenceLooperWarnings( G4bool val);
     // Do not warn (or throw exception) about 'looping' particles
     static G4bool GetSilenceLooperWarnings();
   
  public: // without description    
     static G4bool EnableUseMagneticMoment(G4bool useMoment=true)
     { return EnableMagneticMoment(useMoment); }  // Old name - will be deprecated
   
  public:  // without description

     G4double AtRestGetPhysicalInteractionLength( const G4Track&,
                                                  G4ForceCondition*)
       { return -1.0; }  // No operation in AtRestGPIL

     G4VParticleChange* AtRestDoIt( const G4Track&, const G4Step& )
       { return 0; }     // No operation in AtRestDoIt

     void StartTracking(G4Track* aTrack);
       // Reset state for new (potentially resumed) track 

     virtual void ProcessDescription(std::ostream& outFile) const; // override;
     void PrintStatistics( std::ostream& outStr) const;
   
  protected:

     void SetTouchableInformation(const G4TouchableHandle& touchable);

     void ReportMissingLogger(const char * methodName);
   
  protected:

     G4Navigator* fLinearNavigator;
       // The navigator for the 'mass' geometry
       // (the real one, that physics occurs in)
     G4PropagatorInField* fFieldPropagator;
       // The Propagators used to transport the particle

     G4ThreeVector fTransportEndPosition=     G4ThreeVector( 0.0, 0.0, 0.0 );
     G4ThreeVector fTransportEndMomentumDir=  G4ThreeVector( 0.0, 0.0, 0.0 );
     G4double      fTransportEndKineticEnergy= 0.0;
     G4ThreeVector fTransportEndSpin=  G4ThreeVector( 0.0, 0.0, 0.0 );
     G4bool        fMomentumChanged=   true;
     G4bool        fEndGlobalTimeComputed= false; 
     G4double      fCandidateEndGlobalTime= 0.0;
       // The particle's state after this Step, Store for DoIt

     G4bool        fAnyFieldExists= false; 
   
     G4bool fParticleIsLooping = false;
     G4bool fNewTrack= true;          // Flag from StartTracking 
     G4bool fFirstStepInVolume= true;
     G4bool fLastStepInVolume= false;  // Last step - almost same as next flag
                                // (temporary redundancy for checking) 
     G4bool fGeometryLimitedStep= true;
       // Flag to determine whether a boundary was reached

     G4bool fFieldExertedForce= false; // During current step

     G4TouchableHandle fCurrentTouchableHandle;
     
     G4ThreeVector fPreviousSftOrigin;
     G4double      fPreviousSafety; 
       // Remember last safety origin & value.

     G4ParticleChangeForTransport fParticleChange;
       // New ParticleChange

     G4double fEndPointDistance;

     // Thresholds for looping particles: 
     //
     G4double fThreshold_Warning_Energy =   1.0 * CLHEP::keV;  //  Warn above this energy
     G4double fThreshold_Important_Energy = 1.0 * CLHEP::MeV;  //  Give a few trial above this E
     G4int    fThresholdTrials = 10;       //  Number of trials an important looper survives
       // Above 'important' energy a 'looping' particle in field will 
       // *NOT* be abandoned, except after fThresholdTrials attempts.
     G4int    fAbandonUnstableTrials = 0;  //  Number of trials after which to abandon
                                           //   unstable loopers ( 0 = never )
     // Counter for steps in which particle reports 'looping',
     //  ( Used if it is above 'Important' Energy. )
     G4int    fNoLooperTrials= 0; 

     // Statistics for tracks abandoned due to looping - and 'saved' despite looping
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
     // Whether to avoid calling G4Navigator for short step ( < safety)
     // If using it, the safety estimate for endpoint will likely be smaller.
     //
     G4bool   fShortStepOptimisation; 

     G4SafetyHelper* fpSafetyHelper;    // To pass it the safety value obtained
     G4TransportationLogger* fpLogger;  // Reports issues / raises warnings

  protected:

     static G4bool fUseMagneticMoment;
     static G4bool fUseGravity;
     static G4bool fSilenceLooperWarnings;  // Flag to *Supress* all 'looper' warnings
   
};

#include "G4Transportation.icc"

#endif  
