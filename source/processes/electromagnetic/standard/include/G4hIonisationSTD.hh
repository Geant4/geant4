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
// File name:     G4hIonisationSTD
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
// 07-01-02 new design of em processes (V.Ivanchenko)
// 26-12-02 secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
//
// Class Description:
//
// This class manages the ionisation process for hadrons.
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLoss.
//

// -------------------------------------------------------------------
//

#ifndef G4hIonisationSTD_h
#define G4hIonisationSTD_h 1

#include "G4VEnergyLossSTD.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "globals.hh"

class G4Material;

class G4hIonisationSTD : public G4VEnergyLossSTD
{

public:

  G4hIonisationSTD(const G4String& name = "hIoni");

  ~G4hIonisationSTD();

  G4bool IsApplicable(const G4ParticleDefinition& p)
    {return (p.GetPDGCharge() != 0.0 && p.GetPDGMass() > 10.0*MeV);};

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
                                    const G4Material*, G4double cut);

  virtual G4std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&,
			           G4double&,
			           G4double&,
                                   G4double&);

  virtual void SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double&,
                                   G4double&);

  void SetSubCutoff(G4bool val);

  void PrintInfoDefinition() const;
  // Print out of the class parameters

protected:

  virtual const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle);

private:

  void InitialiseProcess();

  // hide assignment operator
  G4hIonisationSTD & operator=(const G4hIonisationSTD &right);
  G4hIonisationSTD(const G4hIonisationSTD&);

  G4double   mass;
  G4double   ratio;

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;

  G4bool                      subCutoff;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hIonisationSTD::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                   const G4Material*,
                                                         G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double y = electron_mass_c2/mass;
  G4double g = x*y + sqrt((1. + x)*(1. + x*y*y));
  return mass*(g - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hIonisationSTD::MaxSecondaryEnergy(const G4DynamicParticle* dynParticle)
{
  G4double gamma= dynParticle->GetKineticEnergy()/mass + 1.0;
  G4double tmax = 2.0*electron_mass_c2*(gamma*gamma - 1.) /
                  (1. + 2.0*gamma*ratio + ratio*ratio);

  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VSubCutoffProcessor.hh"

inline G4std::vector<G4Track*>*  G4hIonisationSTD::SecondariesAlongStep(
                           const G4Step&   step,
	             	         G4double& tmax,
			         G4double& eloss,
                                 G4double& kinEnergy)
{
  G4std::vector<G4Track*>* newp = 0;
  if(subCutoff) {
    G4VSubCutoffProcessor* sp = SubCutoffProcessor(CurrentMaterialCutsCoupleIndex());
    if (sp) {
      G4VEmModel* model = SelectModel(kinEnergy);
      newp = sp->SampleSecondaries(step,tmax,eloss,model);
    }
  }
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmModel.hh"

inline void G4hIonisationSTD::SecondariesPostStep(
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
  finalP  = finalP.unit();
  aParticleChange.SetMomentumDirectionChange(finalP);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
