// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ionIonisation.hh,v 1.1 1999-01-07 16:11:19 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- G4ionIonisation physics process -----------
//                by Laszlo Urban, 08 Dec 1998 
// ************************************************************
// It is the first implementation of the ionisation for IONS   
// ------------------------------------------------------------
 
#ifndef G4ionIonisation_h
#define G4ionIonisation_h 1
 
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4hEnergyLoss.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4PhysicsLogVector.hh"
#include "G4PhysicsLinearVector.hh"
 
 
class G4ionIonisation : public G4hEnergyLoss 
 
{
  public:
 
     G4ionIonisation(const G4String& processName = "ionIoni"); 

    ~G4ionIonisation();

    G4bool IsApplicable(const G4ParticleDefinition&);

    void PrintInfoDefinition();

    G4double GetMeanFreePath(
                             const G4Track& track,
                             G4double previousStepSize,
                             G4ForceCondition* condition ) ;
 
    G4VParticleChange *PostStepDoIt(const G4Track& track,
                                    const G4Step& Step  ) ;                 

  protected:

    G4double ComputeMicroscopicCrossSection(
                            const G4DynamicParticle* aParticle,
                            G4double KineticEnergy,
                            G4double AtomicNumber);

  private:

  // hide assignment operator 
    G4ionIonisation & operator=(const G4ionIonisation &right);
    G4ionIonisation(const G4ionIonisation&);

  private:
  //  private data members ...............................

    G4double* DeltaCutInKineticEnergy ; 
    G4double  DeltaCutInKineticEnergyNow ;
};
 
#include "G4ionIonisation.icc"
 
#endif
 


