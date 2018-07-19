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
//	G4IonCoulombScatteringModel.hh
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:    G4IonCoulombScatteringModel 
//
// Author:	Cristina Consolandi
//
// Creation date: 05.10.2010 from G4eCoulombScatteringModel 
//				& G4CoulombScatteringModel
//
//
// Class Description:
//	Single Scattering Model for
//      for protons, alpha and heavy Ions
//
// Reference:		
//	M.J. Boschini et al. "Nuclear and Non-Ionizing Energy-Loss 
//	for Coulomb ScatteredParticles from Low Energy up to Relativistic 
//	Regime in Space	Radiation Environment"
//	Accepted for publication in the Proceedings of  the  ICATPP Conference
//	on Cosmic Rays for Particle and Astroparticle Physics, Villa  Olmo, 7-8
//	October,  2010, to be published by World Scientific (Singapore).
//
//      Available for downloading at:
//      http://arxiv.org/abs/1011.4822
//
// -------------------------------------------------------------------
//

#ifndef G4IonCoulombScatteringModel_h
#define G4IonCoulombScatteringModel_h 1

#include "G4VEmModel.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4IonCoulombCrossSection.hh"

#include <vector>

class G4ParticleChangeForGamma;
class G4ParticleDefinition;
class G4IonTable;

class G4IonCoulombScatteringModel : public G4VEmModel
{
public:

  explicit G4IonCoulombScatteringModel(const G4String& nam = 
				       "IonCoulombScattering");
 
  virtual ~G4IonCoulombScatteringModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) final;
 
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


  	
  inline void  SetRecoilThreshold(G4double eth);

  inline void  SetHeavyIonCorr(G4int b);

  inline G4int GetHeavyIonCorr();

  //protected: 
	 
private:

  inline void DefineMaterial(const G4MaterialCutsCouple*);
  
  inline void SetupParticle(const G4ParticleDefinition*);

  // hide assignment operator
  G4IonCoulombScatteringModel & operator=
  (const G4IonCoulombScatteringModel &right) = delete;
  G4IonCoulombScatteringModel(const  G4IonCoulombScatteringModel&) = delete;

  //protected:

  G4IonTable*               theIonTable;
  G4ParticleChangeForGamma* fParticleChange; 
  G4NistManager*            fNistManager;
  G4IonCoulombCrossSection* ioncross;	  

  const std::vector<G4double>* pCuts;
  const G4MaterialCutsCouple* currentCouple; 
  const G4Material*           currentMaterial;
  const G4Element*            currentElement;
  G4int                       currentMaterialIndex;
  
  G4int 		      heavycorr;

  G4double                  cosThetaMin;
  G4double                  recoilThreshold;
				
  // projectile
  const G4ParticleDefinition* particle;		
  const G4ParticleDefinition* theProton;	
  G4double                  mass;		

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void 
G4IonCoulombScatteringModel::DefineMaterial(const G4MaterialCutsCouple* cup)
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
void G4IonCoulombScatteringModel::SetupParticle(const G4ParticleDefinition* p)
{
  if(p != particle) {
    particle = p;
    mass = particle->GetPDGMass();
    ioncross->SetupParticle(p);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4IonCoulombScatteringModel::SetRecoilThreshold(G4double eth)
{
  recoilThreshold = eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4IonCoulombScatteringModel::SetHeavyIonCorr(G4int b) 
{
  heavycorr = b; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4int G4IonCoulombScatteringModel::GetHeavyIonCorr() 
{
  return heavycorr; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
