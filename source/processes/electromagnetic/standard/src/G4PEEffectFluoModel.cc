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
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4PEEffectFluoModel
//
// Author:        Vladimir Ivanchenko on base of G4PEEffectModel
//
// Creation date: 13.06.2010
//
// Modifications:
//
// Class Description:
// Implementation of the photo-electric effect with deexcitation
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4PEEffectFluoModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "Randomize.hh"
#include "G4Material.hh"
#include "G4DataVector.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "G4SauterGavrilaAngularDistribution.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

G4PEEffectFluoModel::G4PEEffectFluoModel(const G4String& nam)
  : G4VEmModel(nam)
{
  theGamma    = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();
  fminimalEnergy = 1.0*CLHEP::eV;
  SetDeexcitationFlag(true);

  fSandiaCof.resize(4,0.0);

  // default generator
  SetAngularDistribution(new G4SauterGavrilaAngularDistribution());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PEEffectFluoModel::~G4PEEffectFluoModel() = default;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PEEffectFluoModel::Initialise(const G4ParticleDefinition*,
				     const G4DataVector&)
{
  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  fPEBelowKShell = G4EmParameters::Instance()->PhotoeffectBelowKShell();
  if(nullptr == fParticleChange) { 
    fParticleChange = GetParticleChangeForGamma(); 
  }
  std::size_t nmat = G4Material::GetNumberOfMaterials();
  fMatEnergyTh.resize(nmat, 0.0);
  for(std::size_t i=0; i<nmat; ++i) { 
    fMatEnergyTh[i] = (*(G4Material::GetMaterialTable()))[i]
      ->GetSandiaTable()->GetSandiaCofForMaterial(0, 0);
    //G4cout << "G4PEEffectFluoModel::Initialise Eth(eV)= " 
    //	   << fMatEnergyTh[i]/eV << G4endl; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4double 
G4PEEffectFluoModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
						G4double energy,
						G4double Z, G4double,
						G4double, G4double)
{
  // This method may be used only if G4MaterialCutsCouple pointer
  //   has been set properly
  CurrentCouple()->GetMaterial()
    ->GetSandiaTable()->GetSandiaCofPerAtom((G4int)Z, energy, fSandiaCof);

  G4double x1 = 1 / energy;

  return x1 * (fSandiaCof[0] + x1 * (fSandiaCof[1] +
    x1 * (fSandiaCof[2] + x1 * fSandiaCof[3])));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4PEEffectFluoModel::CrossSectionPerVolume(const G4Material* material,
					   const G4ParticleDefinition*,
					   G4double energy,
					   G4double, G4double)
{
  // This method may be used only if G4MaterialCutsCouple pointer
  //   has been set properly
  energy = std::max(energy, fMatEnergyTh[material->GetIndex()]);
  const G4double* SandiaCof = 
    material->GetSandiaTable()->GetSandiaCofForMaterial(energy);

  G4double x1 = 1 / energy;

  return x1 * (SandiaCof[0] + x1 * (SandiaCof[1] +
    x1 * (SandiaCof[2] + x1 * SandiaCof[3])));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4PEEffectFluoModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
				       const G4MaterialCutsCouple* couple,
				       const G4DynamicParticle* aDynamicPhoton,
				       G4double,
				       G4double)
{
  SetCurrentCouple(couple);
  const G4Material* aMaterial = couple->GetMaterial();

  G4double energy = aDynamicPhoton->GetKineticEnergy();

  // select randomly one element constituing the material.
  const G4Element* anElement = SelectRandomAtom(aMaterial,theGamma,energy);
  
  //
  // Photo electron
  //

  // Select atomic shell
  G4int nShells = anElement->GetNbOfAtomicShells();
  G4int i = 0;  
  for(; i<nShells; ++i) {
    /*
    G4cout << "i= " << i << " E(eV)= " << energy/eV 
    	   << " Eb(eV)= " << anElement->GetAtomicShell(i)/eV
           << "  " << anElement->GetName() 
	   << G4endl;
    */
    if(energy >= anElement->GetAtomicShell(i)) { break; }
  }

  G4double edep = energy;

  // photo-electron is not sampled if shell is not found or 
  // the flag of photoeffect is "false" and shell is no K 
  if ( (fPEBelowKShell || 0 == i) && i < nShells ) {

    G4double bindingEnergy = anElement->GetAtomicShell(i);
    edep = bindingEnergy;
    G4double esec = 0.0;

    // sample deexcitation cascade
    //
    if(nullptr != fAtomDeexcitation) {
      G4int index = couple->GetIndex();
      if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
	G4int Z = G4lrint(anElement->GetZ());
	auto as = (G4AtomicShellEnumerator)(i);
	const G4AtomicShell* shell = fAtomDeexcitation->GetAtomicShell(Z, as);
        G4double eshell = shell->BindingEnergy();
        if(eshell > bindingEnergy && eshell <= energy) {
          bindingEnergy = eshell;
          edep = eshell;
	}
	std::size_t nbefore = fvect->size();
	fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
	std::size_t nafter = fvect->size();
	for (std::size_t j=nbefore; j<nafter; ++j) {
	  G4double e = ((*fvect)[j])->GetKineticEnergy();
	  if(esec + e > edep) {
	    // correct energy in order to have energy balance
	    e = edep - esec;
	    ((*fvect)[j])->SetKineticEnergy(e);
	    esec += e;
	    /*
	      G4cout << "### G4PEffectFluoModel Edep(eV)= " << edep/eV 
		     << " Esec(eV)= " << esec/eV 
		     << " E["<< j << "](eV)= " << e/eV
		     << " N= " << nafter
		     << " Z= " << Z << " shell= " << i 
		     << "  Ebind(keV)= " << bindingEnergy/keV 
		     << "  Eshell(keV)= " << shell->BindingEnergy()/keV 
		     << G4endl;
	    */
	    // delete the rest of secondaries (should not happens)
	    for (std::size_t jj=nafter-1; jj>j; --jj) { 
	      delete (*fvect)[jj]; 
	      fvect->pop_back(); 
	    }
	    break;	      
	  }
	  esec += e; 
	}
        edep -= esec;
      }
    }
    // create photo electron
    //
    G4double elecKineEnergy = energy - bindingEnergy;
    if (elecKineEnergy > fminimalEnergy) {
      auto aParticle = new G4DynamicParticle(theElectron,
	GetAngularDistribution()->SampleDirection(aDynamicPhoton, 
						  elecKineEnergy,
						  i, couple->GetMaterial()), 
							   elecKineEnergy);
      fvect->push_back(aParticle);
    } else {
      edep += elecKineEnergy;
      elecKineEnergy = 0.0;
    }
    if(std::abs(energy - elecKineEnergy - esec - edep) > CLHEP::eV) {
      G4cout << "### G4PEffectFluoModel dE(eV)= " 
	     << (energy - elecKineEnergy - esec - edep)/eV 
	     << " shell= " << i 
	     << "  E(keV)= " << energy/keV 
	     << "  Ebind(keV)= " << bindingEnergy/keV 
	     << "  Ee(keV)= " << elecKineEnergy/keV 
	     << "  Esec(keV)= " << esec/keV 
	     << "  Edep(keV)= " << edep/keV 
	     << G4endl;
    }
  }

  // kill primary photon
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  if(edep > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
