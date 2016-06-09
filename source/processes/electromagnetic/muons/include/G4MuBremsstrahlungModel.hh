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
// $Id: G4MuBremsstrahlungModel.hh,v 1.13 2005/08/04 08:19:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4MuBremsstrahlungModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 18.05.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 27-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 10-02-04 Add lowestKinEnergy (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
//

//
// Class Description:
//
// Implementation of energy loss for gamma emission by muons

// -------------------------------------------------------------------
//

#ifndef G4MuBremsstrahlungModel_h
#define G4MuBremsstrahlungModel_h 1

#include "G4VEmModel.hh"

class G4Element;
class G4ParticleChangeForLoss;

class G4MuBremsstrahlungModel : public G4VEmModel
{

public:

  G4MuBremsstrahlungModel(const G4ParticleDefinition* p = 0, const G4String& nam = "MuBrem");

  virtual ~G4MuBremsstrahlungModel();

  void SetParticle(const G4ParticleDefinition*);

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  void SetLowestKineticEnergy(G4double e) {lowestKinEnergy = e;};

  G4double MinEnergyCut(const G4ParticleDefinition*,
                        const G4MaterialCutsCouple*);

  G4double ComputeDEDXPerVolume(
                        const G4Material*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy);

  G4double CrossSectionPerVolume(
			const G4Material*,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double cutEnergy,
                              G4double maxEnergy);

  std::vector<G4DynamicParticle*>* SampleSecondaries(
                                const G4MaterialCutsCouple*,
                                const G4DynamicParticle*,
                                      G4double tmin,
                                      G4double maxEnergy);

protected:

  G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
			      G4double kineticEnergy);

public:

  G4double ComputMuBremLoss(G4double Z, G4double A, G4double tkin, G4double cut);

  G4double ComputeMicroscopicCrossSection(G4double tkin,
                                           G4double Z,
                                           G4double A,
                                           G4double cut);

  G4double ComputeDMicroscopicCrossSection(G4double tkin,
                                           G4double Z,
                                           G4double A,
                                           G4double gammaEnergy);

private:

  G4DataVector* ComputePartialSumSigma(const G4Material* material,
                                             G4double tkin, G4double cut);

  const G4Element* SelectRandomAtom(const G4MaterialCutsCouple* couple) const;

  void MakeSamplingTables();


  // hide assignment operator
  G4MuBremsstrahlungModel & operator=(const  G4MuBremsstrahlungModel &right);
  G4MuBremsstrahlungModel(const  G4MuBremsstrahlungModel&);

  G4ParticleDefinition*       theGamma;
  const G4ParticleDefinition* particle;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double lowestKinEnergy;
  G4double minThreshold;
  G4double mass;

  // tables for sampling
  G4int nzdat,ntdat,NBIN;
  static G4double zdat[5],adat[5],tdat[8];
  G4double ya[1001], proba[5][8][1001];
  G4double cutFixed;

  std::vector<G4DataVector*> partialSumSigma;
  G4bool  samplingTablesAreFilled;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4MuBremsstrahlungModel::MaxSecondaryEnergy(
                                 const G4ParticleDefinition*,
    				       G4double kineticEnergy)
{
  return kineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
