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
// File name:     G4eBremsstrahlungSTD
//
// Author:        Laszlo Urban
//
// Creation date: 24.06.1996
//
// Modifications:
//
// 01-10-96 new type G4OrderedTable;  ComputePartialSumSigma()
// 20-03-97 new energy loss+ionisation+brems scheme, L.Urban
// 01-09-98 new method  PrintInfo()
// 10-02-00 modifications , new e.m. structure, L.Urban
// 07-08-00 new cross section/en.loss parametrisation, LPM flag , L.Urban
// 09-08-01 new methods Store/Retrieve PhysicsTable (mma)
// 19-09-01 come back to previous process name "eBrem"
// 29-10-01 all static functions no more inlined (mma)
// 07-01-02 new design of em processes (V.Ivanchenko)
// 26-12-02 secondary production moved to derived classes (VI)
// 24-01-03 Make models region aware (V.Ivanchenko)
//
//
// Class Description:
//
// This class manages the bremsstrahlung for e-/e+
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLoss.
//

// -------------------------------------------------------------------
//

#ifndef G4eBremsstrahlungSTD_h
#define G4eBremsstrahlungSTD_h 1

#include "G4VEnergyLossSTD.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

class G4Material;

class G4eBremsstrahlungSTD : public G4VEnergyLossSTD
{

public:

  G4eBremsstrahlungSTD(const G4String& name = "eBrem");

  ~G4eBremsstrahlungSTD();
 
  G4bool IsApplicable(const G4ParticleDefinition& p) 
    {return (&p == G4Electron::Electron() || &p == G4Positron::Positron());};

  virtual G4double MinPrimaryEnergy(const G4ParticleDefinition*,
                                    const G4Material*, G4double cut)
  {return cut;};

  virtual G4std::vector<G4Track*>* SecondariesAlongStep(
                             const G4Step&, 
			           G4double&,
			           G4double&,
                                   G4double&) {return 0;};

  virtual void SecondariesPostStep(G4ParticleChange&, 
                                   G4VEmModel*, 
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double&,
                                   G4double&);

  void PrintInfoDefinition() const;
  // Print out of the class parameters

protected:

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dynParticle)
  {return dynParticle->GetKineticEnergy();};

private:

  void InitialiseProcess();

  // hide assignment operator
  G4eBremsstrahlungSTD & operator=(const G4eBremsstrahlungSTD &right);
  G4eBremsstrahlungSTD(const G4eBremsstrahlungSTD&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmModel.hh"

inline void G4eBremsstrahlungSTD::SecondariesPostStep(G4ParticleChange& aParticleChange,
                                                      G4VEmModel* model,
                                                const G4MaterialCutsCouple* couple,
                                                const G4DynamicParticle* dp,
                                                      G4double& tcut,
                                                      G4double& kinEnergy)
{
  G4DynamicParticle* gamma = model->SampleSecondary(couple, dp, tcut, kinEnergy);
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.AddSecondary(gamma);
  kinEnergy -= gamma->GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

