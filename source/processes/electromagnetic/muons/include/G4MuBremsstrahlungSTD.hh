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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4MuBremsstrahlungSTD
//
// Author:        Laszlo Urban
// 
// Creation date: 30.09.1997
//
// Modifications:
//
// 10/02/00 modifications , new e.m. structure, L.Urban
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 29-10-01 all static functions no more inlined (mma)
// 10-05-02 V.Ivanchenko update to new design
// 26-12-02 secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
//
// Class Description:
//
// This class manages the Bremsstrahlung process for muons.
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLossSTD.
//

// -------------------------------------------------------------------
//

#ifndef G4MuBremsstrahlungSTD_h
#define G4MuBremsstrahlungSTD_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4VEnergyLossSTD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuBremsstrahlungSTD : public G4VEnergyLossSTD

{
public:

  G4MuBremsstrahlungSTD(const G4String& processName = "MuBrems");

  ~G4MuBremsstrahlungSTD();

  G4bool IsApplicable(const G4ParticleDefinition& p)
            {return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 10.0*MeV);};

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
                                    const G4Material*, G4double cut);

  virtual G4std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&,
			           G4double&,
			           G4double&,
                                   G4double&) {return 0;};

  virtual void SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double&,
                                   G4double&);

  void PrintInfoDefinition() const;
  // Print out of the class parameters

protected:

  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle);

private:

  void InitialiseProcess();


  G4MuBremsstrahlungSTD & operator=(const G4MuBremsstrahlungSTD &right);
  G4MuBremsstrahlungSTD(const G4MuBremsstrahlungSTD&);

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuBremsstrahlungSTD::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                        const G4Material*, 
                                                              G4double cut)
{
  return cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuBremsstrahlungSTD::MaxSecondaryEnergy(const G4DynamicParticle* dynParticle)
{
  return dynParticle->GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmModel.hh"

inline void G4MuBremsstrahlungSTD::SecondariesPostStep(
                                                       G4VEmModel* model,
                                                 const G4MaterialCutsCouple* couple,
                                                 const G4DynamicParticle* dp,
                                                       G4double& tcut,
                                                       G4double& kinEnergy)
{
  G4DynamicParticle* gamma = model->SampleSecondary(couple, dp, tcut, kinEnergy);
  if(gamma) {
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(gamma);
    kinEnergy -= gamma->GetKineticEnergy();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
