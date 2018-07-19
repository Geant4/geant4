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
//	G4eSingleCoulombScatteringModel.hh
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4eSingleCoulombScatteringModel
//
// Author:	Cristina Consolandi
//
// Creation date: 20.10.2011
//
// Class Description:
//      Single Scattering model for electron-nuclei interaction.
//      Suitable for high energy electrons and low scattering angles.
//
// Reference:
//    	M.J. Boschini et al.
//      "Non Ionizing Energy Loss induced by Electrons in the Space Environment"
//       Proc. of the 13th International Conference on Particle Physics and Advanced Technology
//       (13th ICPPAT, Como 3-7/10/2011), World Scientific (Singapore).
//      Available at: http://arxiv.org/abs/1111.4042v4
//
//
// -------------------------------------------------------------------
//

#ifndef G4eSingleCoulombScatteringModel_h
#define G4eSingleCoulombScatteringModel_h 1

#include "G4VEmModel.hh"
#include "globals.hh"
#include "G4ScreeningMottCrossSection.hh"

#include <vector>

class G4ParticleChangeForGamma;
class G4ParticleDefinition;
class G4IonTable;
class G4Nistmanager;

class G4eSingleCoulombScatteringModel : public G4VEmModel
{

public:

  explicit G4eSingleCoulombScatteringModel(const G4String& nam = "eSingleCoulombScat");

  virtual ~G4eSingleCoulombScatteringModel();

  virtual void Initialise(const G4ParticleDefinition*,
			  const G4DataVector&) final;

  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel* masterModel) final;

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
				G4double kinEnergy,
				G4double Z,
				G4double A,
				G4double cut,
				G4double emax) final;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) final;

  inline void SetRecoilThreshold(G4double eth);
  inline void SetXSectionModel(const G4String& model);

private:

  inline void DefineMaterial(const G4MaterialCutsCouple*);
  inline void SetupParticle(const G4ParticleDefinition*);

  // hide assignment operator
  G4eSingleCoulombScatteringModel & operator=
  (const G4eSingleCoulombScatteringModel &right) = delete;
  G4eSingleCoulombScatteringModel(const  G4eSingleCoulombScatteringModel&) = delete;

  G4IonTable*               theIonTable;
  G4ParticleChangeForGamma* fParticleChange;
  G4NistManager*            fNistManager;
  G4ScreeningMottCrossSection* Mottcross;

  const std::vector<G4double>* pCuts;
  const G4MaterialCutsCouple* currentCouple;
  const G4Material*           currentMaterial;
  const G4Element*            currentElement;
  G4int                       currentMaterialIndex;

  G4double                  cosThetaMin;
  G4double                  recoilThreshold;
  G4int                     FormFactor;
  G4int                     XSectionModel;

  // projectile
  const G4ParticleDefinition* particle;
  G4double                  mass;
  G4double                  lowEnergyLimit;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void
G4eSingleCoulombScatteringModel::DefineMaterial(const G4MaterialCutsCouple* cup)
{
  if(cup != currentCouple) {
    currentCouple = cup;
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void
G4eSingleCoulombScatteringModel::SetupParticle(const G4ParticleDefinition* p)
{
  if(p != particle) {
    particle = p;
    mass = particle->GetPDGMass();
    Mottcross->SetupParticle(p);
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eSingleCoulombScatteringModel::SetRecoilThreshold(G4double eth)
{
  recoilThreshold = eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eSingleCoulombScatteringModel::SetXSectionModel(const G4String& model)
{
  if(model == "fast") { XSectionModel=1; }
  else if(model == "precise") { XSectionModel=0; }
  else { G4cout<<"G4eSingleCoulombScatteringModel WARNING: "<<model
	       <<" : G4eSingleScatteringModel x-section model is not valid"<<G4endl;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
