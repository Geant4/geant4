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
// File name:     G4MuPairProductionModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 18.05.2002
//
// Modifications: 
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
//

//
// Class Description:
//
// Implementation of e+e- pair production by muons
//

// -------------------------------------------------------------------
//

#ifndef G4MuPairProductionModel_h
#define G4MuPairProductionModel_h 1

#include "G4VEmModel.hh"

class G4MuPairProductionModel : public G4VEmModel
{

public:

  G4MuPairProductionModel(const G4ParticleDefinition* p = 0);

  ~G4MuPairProductionModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);
  
  G4double HighEnergyLimit(const G4ParticleDefinition* p);

  G4double LowEnergyLimit(const G4ParticleDefinition* p);

  void SetHighEnergyLimit(G4double e) {highKinEnergy = e;};

  void SetLowEnergyLimit(G4double e) {lowKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);

  G4bool IsInCharge(const G4ParticleDefinition*);

  G4double ComputeDEDX(const G4Material*,
                       const G4ParticleDefinition*,
                             G4double kineticEnergy,
                             G4double cutEnergy);

  G4double CrossSection(const G4Material*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy,
                              G4double maxEnergy);

  G4DynamicParticle* SampleSecondary(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

  G4std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

  virtual G4double MaxSecondaryEnergy(
				const G4DynamicParticle* dynParticle);
protected:

  virtual G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
    				            G4double kineticEnergy);

private:

  G4double ComputMuPairLoss(G4double Z, G4double tkin, G4double cut);

  G4double ComputeMicroscopicCrossSection(G4double tkin,
                                          G4double Z,
                                          G4double A,
                                          G4double cut);

  G4double ComputeDMicroscopicCrossSection(G4double tkin,
                                           G4double Z,
                                           G4double pairEnergy);

  G4double ComputeDDMicroscopicCrossSection(G4double tkin,
                                           G4double Z,
                                           G4double pairEnergy,
                                           G4double asymmetry);

  G4DataVector* ComputePartialSumSigma(const G4Material* material,
                                             G4double tkin, G4double cut);

  const G4Element* SelectRandomAtom(const G4MaterialCutsCouple* couple) const;

  void MakeSamplingTables();

  // hide assignment operator
  G4MuPairProductionModel & operator=(const  G4MuPairProductionModel &right);
  G4MuPairProductionModel(const  G4MuPairProductionModel&);

  G4double minPairEnergy;
  G4double highKinEnergy;
  G4double lowKinEnergy;

  // tables for sampling
  G4int nzdat,ntdat,NBIN;
  static G4double zdat[5],adat[5],tdat[8];
  G4double ya[1001],proba[5][8][1001];

  G4std::vector<G4DataVector*> partialSumSigma;
  G4bool  samplingTablesAreFilled;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4MuPairProductionModel::MaxSecondaryEnergy(
				 const G4DynamicParticle* dynParticle)
{
  return dynParticle->GetKineticEnergy();
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4MuPairProductionModel::MaxSecondaryEnergy(
                                 const G4ParticleDefinition*,
    				       G4double kineticEnergy)
{
  return kineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
