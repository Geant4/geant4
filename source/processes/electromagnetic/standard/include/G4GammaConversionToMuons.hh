//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4GammaConversionToMuons.hh,v 1.1 2002-04-19 14:44:32 hbu Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//         ------------ G4GammaConversionToMuons physics process ------
//         by H.Burkhardt, S. Kelner and R. Kokoulin, April 2002
// -----------------------------------------------------------------------------

// class description
//
// gamma ---> mu+ mu-
// inherit from G4VDiscreteProcess
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4GammaConversionToMuons_h
#define G4GammaConversionToMuons_h 1

#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4VDiscreteProcess.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsLogVector.hh"
#include "G4Element.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Step.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4GammaConversionToMuons : public G4VDiscreteProcess

{
  public:  // with description

     G4GammaConversionToMuons(const G4String& processName ="GammaToMuPair");

    ~G4GammaConversionToMuons();

     G4bool IsApplicable(const G4ParticleDefinition&);
       // true for Gamma only.

     void BuildPhysicsTable(const G4ParticleDefinition&);
       // here dummy, the total cross section parametrization is used rather
       // than tables,  just calling PrintInfoDefinition

     void PrintInfoDefinition();
       // Print few lines of informations about the process: validity range,
       // origine ..etc..
       // Invoked by BuildThePhysicsTable().

     void SetCrossSecFactor(G4double fac);
       // Set the factor to artificially increase the cross section (by default 1)

     G4double inline GetCrossSecFactor() { return CrossSecFactor;}
       // Get the factor to artificially increase the cross section

     G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* condition);
       // It returns the MeanFreePath of the process for the current track :
       // (energy, material)
       // The previousStepSize and G4ForceCondition* are not used.
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.

     G4double GetCrossSectionPerAtom(const G4DynamicParticle* aDynamicGamma,
                                           G4Element*         anElement);
       // It returns the total CrossSectionPerAtom of the process,
       // for the current DynamicGamma (energy), in anElement.

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
                                    const G4Step& aStep);
       // It computes the final state of the process (at end of step),
       // returned as a ParticleChange object.
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.

  protected:

     virtual G4double ComputeCrossSectionPerAtom(G4double GammaEnergy,
               G4double AtomicZ,G4double AtomicA);

     virtual G4double ComputeMeanFreePath (G4double GammaEnergy,
                                           G4Material* aMaterial);

  private:

     G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicGamma,
                                 G4Material* aMaterial);

  private:

     // hide assignment operator as private
     G4GammaConversionToMuons& operator=(const G4GammaConversionToMuons &right);
     G4GammaConversionToMuons(const G4GammaConversionToMuons& );

  private:

     G4double LowestEnergyLimit ;     // low  energy limit of the tables
     G4double HighestEnergyLimit ;    // high energy limit of the tables

     G4double fminimalEnergy;         // minimalEnergy of produced particles

     G4double MeanFreePath;           // actual MeanFreePath (current medium)
     G4double CrossSecFactor;         // factor to artificially increase the cross section

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4GammaConversionToMuons.icc"

#endif

