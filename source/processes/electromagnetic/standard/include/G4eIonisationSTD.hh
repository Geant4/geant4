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
// File name:     G4eIonisationSTD
//
// Author:        Laszlo Urban
// 
// Creation date: 20.03.1997
//
// Modifications: 
//
// 10-02-00 modifications , new e.m. structure, L.Urban
// 03-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 13-08-01 new function ComputeRestrictedMeandEdx() (mma)
// 19-09-01 come back to previous ProcessName "eIoni"
// 29-10-01 all static functions no more inlined (mma)
// 07-01-02 new design of em processes (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (VI)
//
//
// Class Description: 
//
// This class manages the ionisation process for e-/e+
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLoss.
// 

// -------------------------------------------------------------------
//

#ifndef G4eIonisationSTD_h
#define G4eIonisationSTD_h 1

#include "G4VEnergyLossSTD.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4Material;
class G4ParticleDefinition;

class G4eIonisationSTD : public G4VEnergyLossSTD
{

public:

  G4eIonisationSTD(const G4String& name = "eIoni");

  ~G4eIonisationSTD();
 
  G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut);

  virtual G4std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&, 
                             const G4Material*, 
                             const G4DynamicParticle*,
			           G4double,
                                   G4double);

  virtual void SecondariesPostStep(G4ParticleChange&, 
                                   G4VEmModel*, 
                             const G4Material*, 
                             const G4DynamicParticle*,
                                   G4double&,
                                   G4double&);

  void SetSubCutoffProcessor(G4VSubCutoffProcessor*);

  G4VSubCutoffProcessor* SubCutoffProcessor() {return subCutoffProcessor;};

  void PrintInfoDefinition() const;
  // Print out of the class parameters

protected:

  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dp);

private:

  void InitialiseProcess();

  // hide assignment operator 
  G4eIonisationSTD & operator=(const G4eIonisationSTD &right);
  G4eIonisationSTD(const G4eIonisationSTD&);

  const G4ParticleDefinition* theElectron;
  G4bool isElectron;  
  G4VSubCutoffProcessor* subCutoffProcessor;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eIonisationSTD::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                   const G4Material*, 
                                                         G4double cut)
{
  G4double x = cut;
  if(isElectron) x += cut;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4eIonisationSTD::IsApplicable(const G4ParticleDefinition& p) 
{
  return (&p == G4Electron::Electron() || &p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eIonisationSTD::MaxSecondaryEnergy(const G4DynamicParticle* dp) 
{
  G4double tmax = dp->GetKineticEnergy();
  if(isElectron) tmax *= 0.5;
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VSubCutoffProcessor.hh"

inline G4std::vector<G4Track*>*  G4eIonisationSTD::SecondariesAlongStep(
                             const G4Step& step, 
                             const G4Material*, 
                             const G4DynamicParticle* dp,
			           G4double tmax,
                                   G4double eloss)
{
  G4std::vector<G4Track*>* newp = 0;
  if(subCutoffProcessor) {
    newp = subCutoffProcessor->SampleSecondaries(step,dp,tmax,eloss);
  }
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

#include "G4VEmModel.hh"

inline void G4eIonisationSTD::SecondariesPostStep(G4ParticleChange& aParticleChange, 
                                                  G4VEmModel* model, 
                                            const G4Material* material, 
                                            const G4DynamicParticle* dp,
                                                  G4double& tcut,
                                                  G4double& kinEnergy)
{
  G4DynamicParticle* delta = model->SampleSecondary(material, dp, tcut, kinEnergy);
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(delta);
  G4ThreeVector finalP = dp->GetMomentum();
  kinEnergy -= delta->GetKineticEnergy();
  finalP -= delta->GetMomentum();
  finalP = finalP.unit();
  aParticleChange.SetMomentumDirectionChange(finalP);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
