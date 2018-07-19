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
// $Id: G4AnnihiToMuPair.hh 97391 2016-06-02 10:08:45Z gcosmo $
//
//         ------------ G4AnnihiToMuPair physics process ------
//         by H.Burkhardt, S. Kelner and R. Kokoulin, November 2002
// -----------------------------------------------------------------------------

// class description
//
// (high energy) e+ (atomic) e- ---> mu+ mu-
// inherit from G4VDiscreteProcess
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......//
//
// 04.02.03 : cosmetic simplifications (mma)
// 27.01.03 : first implementation (hbu)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4AnnihiToMuPair_h
#define G4AnnihiToMuPair_h 1

#include "G4VDiscreteProcess.hh"
#include "globals.hh"

class G4ParticleDefinition;
class G4Track;
class G4Step;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4AnnihiToMuPair : public G4VDiscreteProcess
{
  public:  // with description

     explicit G4AnnihiToMuPair(const G4String& processName ="AnnihiToMuPair",
		            G4ProcessType type = fElectromagnetic);

    ~G4AnnihiToMuPair();

     G4bool IsApplicable(const G4ParticleDefinition&) override;
       // true for positron only.

     void BuildPhysicsTable(const G4ParticleDefinition&) override;
       // here dummy, just calling PrintInfoDefinition
       // the total cross section is calculated analytically

     void PrintInfoDefinition();
       // Print few lines of informations about the process: validity range,
       // origine ..etc..
       // Invoked by BuildPhysicsTable().

     void SetCrossSecFactor(G4double fac);
       // Set the factor to artificially increase the crossSection (default 1)

     G4double GetCrossSecFactor() {return CrossSecFactor;};
       // Get the factor to artificially increase the cross section

     G4double CrossSectionPerVolume(G4double PositronEnergy, 
				    const G4Material*);
       // Compute total cross section					 

     G4double ComputeCrossSectionPerAtom(G4double PositronEnergy,
                                         G4double AtomicZ);
       // Compute total cross section					 

     G4double GetMeanFreePath(const G4Track& aTrack,
                              G4double previousStepSize,
                              G4ForceCondition* ) override;
       // It returns the MeanFreePath of the process for the current track :
       // (energy, material)
       // The previousStepSize and G4ForceCondition* are not used.
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.

     G4VParticleChange* PostStepDoIt(const G4Track& aTrack,
				     const G4Step& aStep) override;
       // It computes the final state of the process (at end of step),
       // returned as a ParticleChange object.
       // This function overloads a virtual function of the base class.
       // It is invoked by the ProcessManager of the Particle.

  private:

     // hide assignment operator as private
     G4AnnihiToMuPair& operator=(const G4AnnihiToMuPair &right) = delete;
     G4AnnihiToMuPair(const G4AnnihiToMuPair& ) = delete;

     G4double LowestEnergyLimit;     // Energy threshold of e+
     G4double HighestEnergyLimit;    // Limit of validity of the model

     G4double CurrentSigma;          // the last value of cross section per volume 

     G4double CrossSecFactor;        // factor to artificially increase 
                                     // the cross section, static to make sure
				     // to have single value
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

