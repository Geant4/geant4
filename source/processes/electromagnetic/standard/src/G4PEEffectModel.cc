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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PEEffectModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 21.03.2005
//
// Modifications:
//
// 04.12.05 : SetProposedKineticEnergy(0.) for the killed photon (mma)
// 20.02.09 : Added initialisation of deexcitation flag and method
//            CrossSectionPerVolume instead of mfp (V.Ivanchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4PEEffectModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4SauterGavrilaAngularDistribution.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4PEEffectModel::G4PEEffectModel(const G4ParticleDefinition*, 
				 const G4String& nam)
  : G4VEmModel(nam)
{
  G4cout << "### G4PEEffectModel is obsolete "
	 << "and will be removed for the next release." << G4endl;

  theGamma    = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  fminimalEnergy = 1.0*eV;
  fParticleChange = 0;

  // default generator
  SetAngularDistribution(new G4SauterGavrilaAngularDistribution());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PEEffectModel::~G4PEEffectModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PEEffectModel::Initialise(const G4ParticleDefinition*,
				 const G4DataVector&)
{
  // always false before the run
  SetDeexcitationFlag(false);
  if(!fParticleChange) { fParticleChange = GetParticleChangeForGamma(); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double G4PEEffectModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
						     G4double energy,
						     G4double Z, G4double,
						     G4double, G4double)
{
  G4double* SandiaCof = G4SandiaTable::GetSandiaCofPerAtom((G4int)Z, energy);

  G4double energy2 = energy*energy;
  G4double energy3 = energy*energy2;
  G4double energy4 = energy2*energy2;

  return SandiaCof[0]/energy  + SandiaCof[1]/energy2 +
    SandiaCof[2]/energy3 + SandiaCof[3]/energy4;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PEEffectModel::CrossSectionPerVolume(const G4Material* material,
						const G4ParticleDefinition*,
						G4double energy,
						G4double, G4double)
{
  G4double* SandiaCof = 
    material->GetSandiaTable()->GetSandiaCofForMaterial(energy);
				
  G4double energy2 = energy*energy;
  G4double energy3 = energy*energy2;
  G4double energy4 = energy2*energy2;
	  
  return SandiaCof[0]/energy  + SandiaCof[1]/energy2 +
    SandiaCof[2]/energy3 + SandiaCof[3]/energy4; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PEEffectModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					const G4MaterialCutsCouple* couple,
					const G4DynamicParticle* aDynamicPhoton,
					G4double,
					G4double)
{
  const G4Material* aMaterial = couple->GetMaterial();

  G4double energy = aDynamicPhoton->GetKineticEnergy();
  G4ParticleMomentum PhotonDirection = aDynamicPhoton->GetMomentumDirection();

  // select randomly one element constituing the material.
  const G4Element* anElement = SelectRandomAtom(aMaterial,theGamma,energy);
  
  //
  // Photo electron
  //

  // Select atomic shell
  G4int nShells = anElement->GetNbOfAtomicShells();
  G4int i  = 0;  
  while ((i<nShells) && (energy<anElement->GetAtomicShell(i))) { ++i; }

  G4double edep = energy;

  // no shell available
  if (i < nShells) {
  
    G4double bindingEnergy = anElement->GetAtomicShell(i);
    G4double elecKineEnergy = energy - bindingEnergy;

    if (elecKineEnergy > fminimalEnergy) {

      edep = bindingEnergy;
      G4ThreeVector elecDirection =
	GetAngularDistribution()->SampleDirection(aDynamicPhoton, 
						  elecKineEnergy,
						  i, 
						  couple->GetMaterial());

      G4DynamicParticle* aParticle = 
	new G4DynamicParticle(theElectron, elecDirection, elecKineEnergy);
      fvect->push_back(aParticle);

    } 
  }
    
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  if(edep > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
