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
// File name:     G4MuIonisationSTD
//
// Author:        Laszlo Urban
// 
// Creation date: 30.05.1997
//
// Modifications: 
//
// corrected by L.Urban on 24/09/97
// corrected by L.Urban on 13/01/98
// bugs fixed by L.Urban on 02/02/99
// 10/02/00 modifications , new e.m. structure, L.Urban
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 14-08-01 new function ComputeRestrictedMeandEdx() + 'cleanup' (mma)
// 19-09-01 come back to previous process name "hIoni"
// 29-10-01 all static functions no more inlined
// 10-05-02 V.Ivanchenko update to new design
// 09-12-02 V.Ivanchenko remove warning
// 26-12-02 Secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
//
// Class Description:
//
// This class manages the ionisation process for muons.
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLossSTD.
//

// -------------------------------------------------------------------
//

#ifndef G4MuIonisationSTD_h
#define G4MuIonisationSTD_h 1

#include "G4VEnergyLossSTD.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"

class G4Material;

class G4MuIonisationSTD : public G4VEnergyLossSTD
{

public:

  G4MuIonisationSTD(const G4String& name = "muIoni");

  ~G4MuIonisationSTD();

  G4bool IsApplicable(const G4ParticleDefinition& p)
    {return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 10.0*MeV);};

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
                                    const G4Material*, G4double cut);

  virtual G4std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&,
			           G4double&,
			           G4double&,
                                   G4double&);

  virtual void SecondariesPostStep(G4ParticleChange&,
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double&,
                                   G4double&);

  void SetSubCutoffProcessor(G4VSubCutoffProcessor*);

  G4VSubCutoffProcessor* SubCutoffProcessor() {return subCutoffProcessor;};

  void SetSubCutoff(G4bool val);

  void PrintInfoDefinition() const;
  // Print out of the class parameters

protected:

  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle);

private:

  void InitialiseProcess();

  // hide assignment operator 
  G4MuIonisationSTD & operator=(const G4MuIonisationSTD &right);
  G4MuIonisationSTD(const G4MuIonisationSTD&);

  G4double    mass;
  G4double    ratio;

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;

  G4VSubCutoffProcessor* subCutoffProcessor;
  G4bool                 subCutoff;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuIonisationSTD::MinPrimaryEnergy(const G4ParticleDefinition* p,
                                                    const G4Material*, 
                                                          G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double y = electron_mass_c2/mass;
  G4double g = x*y + sqrt((1. + x)*(1. + x*y*y));
  return mass*(g - 1.0); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuIonisationSTD::MaxSecondaryEnergy(const G4DynamicParticle* dynParticle)
{
  G4double gamma= dynParticle->GetKineticEnergy()/mass + 1.0;
  G4double tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);
  
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VSubCutoffProcessor.hh"

inline G4std::vector<G4Track*>*  G4MuIonisationSTD::SecondariesAlongStep(
                           const G4Step&   step, 
	             	         G4double& tmax,
			         G4double& eloss,
                                 G4double& kinEnergy)
{
  G4std::vector<G4Track*>* newp = 0;
  if(subCutoff) {
    G4VEmModel* model = SelectModel(kinEnergy);
    newp = subCutoffProcessor->SampleSecondaries(step,tmax,eloss,model);
  }
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmModel.hh"

inline void G4MuIonisationSTD::SecondariesPostStep(G4ParticleChange& aParticleChange,
                                                   G4VEmModel* model,
                                             const G4MaterialCutsCouple* couple,
                                             const G4DynamicParticle* dp,
                                                   G4double& tcut,
                                                   G4double& kinEnergy)
{
  G4DynamicParticle* delta = model->SampleSecondary(couple, dp, tcut, kinEnergy);
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(delta);
  G4ThreeVector finalP = dp->GetMomentum();
  kinEnergy -= delta->GetKineticEnergy();
  finalP -= delta->GetMomentum();
  finalP = finalP.unit();
  aParticleChange.SetMomentumDirectionChange(finalP);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
