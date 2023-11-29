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

#include "G4Transportation.hh"

#include "G4Track.hh"
#include "G4Step.hh"

class G4PathFinder;

class G4CoupledTransportation : public G4Transportation
{

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

     // AlongStepDoIt is implemented by G4Transportation.

     G4VParticleChange* PostStepDoIt(
                             const G4Track& track,
                             const G4Step&  stepData
                            );
       // Responsible for the relocation

     // PostStepGetPhysicalInteractionLength is implemented by
     // G4Transportation to force PostStepDoIt, but not limiting the step.

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
     G4bool IsFirstStepInAnyVolume() const { return fFirstStepInVolume; }
     G4bool IsLastStepInAnyVolume() const { return fGeometryLimitedStep; }
     G4bool IsFirstStepInMassVolume() const { return fFirstStepInMassVolume; }
     G4bool IsLastStepInMassVolume() const { return fMassGeometryLimitedStep; } 

  public:  // without description

     void StartTracking(G4Track* aTrack); 
     void EndTracking();

     static G4bool EnableUseMagneticMoment(G4bool useMoment=true)
      { return EnableMagneticMoment(useMoment); }
     // Old name ... obsolete
   
  protected:

     void ReportInexactEnergy(G4double startEnergy, G4double endEnergy);
       // Issue warning

     void ReportMove( G4ThreeVector OldVector, G4ThreeVector NewVector,
                      const G4String& Quantity );
   
  private:

     G4PathFinder*        fPathFinder;
       // The PathFinder used to transport the particle

     G4double             fPreviousMassSafety;
     G4double             fPreviousFullSafety;

     G4bool fMassGeometryLimitedStep;
       // Flag to determine whether a 'mass' boundary was reached.

  private:

     G4bool fFirstStepInMassVolume;
     // G4bool fLastStepInMassVolume; => use fMassGeometryLimitedStep 

     static G4bool fSignifyStepInAnyVolume;
       // True: First/Last step in any one of the geometries
       // False: First/Last step in volume of 'mass' geometry
};

#endif  
