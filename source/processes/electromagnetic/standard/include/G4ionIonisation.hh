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
// $Id: G4ionIonisation.hh,v 1.33 2004/12/01 19:37:13 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4ionIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 07.05.2002
//
// Modifications:
//
// 26-12-02 Secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 13-02-03 SubCutoff regime is assigned to a region (V.Ivanchenko)
// 15-02-03 Add control on delta pointer (V.Ivanchenko)
// 23-05-03 Add fluctuation model as a member function (V.Ivanchenko)
// 03-08-03 Add effective charge and saturation of tmax (V.Ivanchenko)
// 12-11-03 Fix problem of negative effective charge (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 21-01-04 Migrade to G4ParticleChangeForLoss (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
//
// Class Description:
//
// This class manages the ionisation process for ions.
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLoss.
//

// -------------------------------------------------------------------
//

#ifndef G4ionIonisation_h
#define G4ionIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "G4ionEffectiveCharge.hh"

class G4Material;
class G4VEmFluctuationModel;

class G4ionIonisation : public G4VEnergyLossProcess
{

public:

  G4ionIonisation(const G4String& name = "ionIoni");

  virtual ~G4ionIonisation();

  G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition* p,
                                    const G4Material*, G4double cut);

  virtual std::vector<G4Track*>* SecondariesAlongStep(
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

  void PrintInfoDefinition();
  // Print out of the class parameters

protected:

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
                                           const G4ParticleDefinition*);

  virtual G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle);

private:

  // hide assignment operator
  G4ionIonisation & operator=(const G4ionIonisation &right);
  G4ionIonisation(const G4ionIonisation&);

  G4ionEffectiveCharge        effCharge;

  // cash
  const G4Material*           theMaterial;
  const G4ParticleDefinition* currentParticle;
  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;
  G4VEmFluctuationModel*      flucModel;

  G4bool                      isInitialised;
  G4bool                      subCutoff;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4ionIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() != 0.0 && !p.IsShortLived() && p.GetParticleType() == "nucleus");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ionIonisation::MinPrimaryEnergy(
          const G4ParticleDefinition*, const G4Material*, G4double cut)
{
  G4double x = 0.5*cut/electron_mass_c2;
  G4double g = std::sqrt(1. + x);
  return proton_mass_c2*(g - 1.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4ionIonisation::MaxSecondaryEnergy(const G4DynamicParticle* dynParticle)
{
  G4double mass  = dynParticle->GetMass();
  G4double gamma = dynParticle->GetKineticEnergy()/mass + 1.0;
  G4double tmax  = electron_mass_c2*std::min(2.0*(gamma*gamma - 1.),
                                             51200.*std::pow(proton_mass_c2/mass,0.66667));
  return tmax;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


inline G4double G4ionIonisation::GetMeanFreePath(const G4Track& track,
                                                       G4double, 
                                                       G4ForceCondition* cond)
{

  currentParticle = track.GetDefinition();
  theMaterial     = track.GetMaterial();
  G4double q_2    = effCharge.EffectiveChargeSquareRatio(currentParticle,theMaterial,
                                                         track.GetKineticEnergy());
  G4double mRatio = proton_mass_c2/track.GetDynamicParticle()->GetMass();

  SetMassRatio(mRatio);
  SetReduceFactor(1.0/(q_2*mRatio));
  SetChargeSquare(q_2);
  SetChargeSquareRatio(q_2);

  return G4VEnergyLossProcess::GetMeanFreePath(track, 0.0, cond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VSubCutoffProcessor.hh"

inline std::vector<G4Track*>*  G4ionIonisation::SecondariesAlongStep(
                           const G4Step&   step,
	             	         G4double& tmax,
			         G4double& eloss,
                                 G4double& scaledEnergy)
{
  std::vector<G4Track*>* newp = 0;
  G4double e = step.GetTrack()->GetKineticEnergy() - eloss;
  fParticleChange.SetProposedCharge(effCharge.EffectiveCharge(currentParticle,theMaterial,e));
  if(subCutoff) {
    G4VSubCutoffProcessor* sp = SubCutoffProcessor(CurrentMaterialCutsCoupleIndex());
    if (sp) {
      G4VEmModel* model = SelectModel(scaledEnergy);
      newp = sp->SampleSecondaries(step,tmax,eloss,model);
    }
  }
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmModel.hh"

inline void G4ionIonisation::SecondariesPostStep(
                                                 G4VEmModel* model,
                                           const G4MaterialCutsCouple* couple,
                                           const G4DynamicParticle* dp,
                                                 G4double& tcut,
                                                 G4double& kinEnergy)
{
  G4DynamicParticle* delta = model->SampleSecondary(couple, dp, tcut, kinEnergy);
  if(delta) {
    fParticleChange.SetNumberOfSecondaries(1);
    fParticleChange.AddSecondary(delta);
    G4ThreeVector finalP = dp->GetMomentum();
    kinEnergy -= delta->GetKineticEnergy();
    fParticleChange.SetProposedCharge(effCharge.EffectiveCharge(currentParticle,theMaterial,kinEnergy));
    finalP -= delta->GetMomentum();
    finalP = finalP.unit();
    fParticleChange.SetProposedMomentumDirection(finalP);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
