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
// File name:     G4MuPairProductionSTD
//
// Author:        Laszlo Urban
// 
// Creation date: 30.05.1998
//
// Modifications: 

// 10-02-00 modifications , new e.m. structure, L.Urban
// 10-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 29-10-01 all static functions no more inlined (mma)
// 10-05-02 V.Ivanchenko update to new design
// 26-12-02 secondary production moved to derived classes (VI)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
//
// Class Description:
//
// This class manages the PairProduction process for muons.
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLossSTD.
//

// -------------------------------------------------------------------
//

#ifndef G4MuPairProductionSTD_h
#define G4MuPairProductionSTD_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4VEnergyLossSTD.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuPairProductionSTD : public G4VEnergyLossSTD

{
public:

  G4MuPairProductionSTD(const G4String& processName = "MuPairProd");

  ~G4MuPairProductionSTD();

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

  void SetSubCutoffProcessor(G4VSubCutoffProcessor*);

  G4VSubCutoffProcessor* SubCutoffProcessor() {return subCutoffProcessor;};

  void SetSubCutoff(G4bool val);

  void PrintInfoDefinition() const;
  // Print out of the class parameters

protected:

  const G4ParticleDefinition* DefineBaseParticle(const G4ParticleDefinition* p);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dp);

private:

  void InitialiseProcess();


  G4MuPairProductionSTD & operator=(const G4MuPairProductionSTD &right);
  G4MuPairProductionSTD(const G4MuPairProductionSTD&);

  const G4ParticleDefinition* theParticle;
  const G4ParticleDefinition* theBaseParticle;

  G4VSubCutoffProcessor* subCutoffProcessor;
  G4bool                      subCutoff;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuPairProductionSTD::MinPrimaryEnergy(const G4ParticleDefinition*,
                                                        const G4Material*,
                                                              G4double cut)
{
  return G4std::max(cut, 2.0*electron_mass_c2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuPairProductionSTD::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{
  return dp->GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VSubCutoffProcessor.hh"

inline G4std::vector<G4Track*>*  G4MuPairProductionSTD::SecondariesAlongStep(
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

inline void G4MuPairProductionSTD::SecondariesPostStep(
                                   G4VEmModel* model,
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp,
                                   G4double& tcut,
                                   G4double& kinEnergy)
{
  G4std::vector<G4DynamicParticle*>* newp =
         model->SampleSecondaries(couple, dp, tcut, kinEnergy);
  if(newp) {
    aParticleChange.SetNumberOfSecondaries(2);
    G4DynamicParticle* elpos = (*newp)[0];
    aParticleChange.AddSecondary(elpos);
    kinEnergy -= elpos->GetKineticEnergy();
    elpos = (*newp)[1];
    aParticleChange.AddSecondary(elpos);
    kinEnergy -= elpos->GetKineticEnergy();
    delete newp;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
