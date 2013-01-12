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
// $Id: G4LivermorePolarizedPhotoElectricModel.hh,v 1.2 2010-11-23 16:42:15 flongo Exp $
//
// Authors: G.Depaola & F.Longo
//

#include "G4LivermorePolarizedPhotoElectricModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"
#include "G4SauterGavrilaAngularDistribution.hh"
#include "G4AtomicShell.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricModel::G4LivermorePolarizedPhotoElectricModel(
   const G4String& nam)
  :G4VEmModel(nam),fParticleChange(0),
   crossSectionHandler(0), shellCrossSectionHandler(0)
{
  lowEnergyLimit  = 10 * eV; // SI - Could be 10 eV ?
  highEnergyLimit = 100 * GeV;
  
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  theGamma    = G4Gamma::Gamma();
  theElectron = G4Electron::Electron();

  // use default generator
  SetAngularDistribution(new G4SauterGavrilaAngularDistribution());

  // polarized generator needs fix
  //SetAngularDistribution(new G4PhotoElectricAngularGeneratorPolarized());
  SetDeexcitationFlag(true);
  fAtomDeexcitation = 0; 
  fDeexcitationActive = false; 

  if (verboseLevel > 1) {
    G4cout << "Livermore Polarized PhotoElectric is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / keV << " keV - "
	   << highEnergyLimit / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePolarizedPhotoElectricModel::~G4LivermorePolarizedPhotoElectricModel()
{  
  delete crossSectionHandler;
  delete shellCrossSectionHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePolarizedPhotoElectricModel::Initialise(
                    const G4ParticleDefinition*,
		    const G4DataVector&)
{
  if (verboseLevel > 3) {
    G4cout << "Calling G4LivermorePolarizedPhotoElectricModel::Initialise()" << G4endl;
  }
  if (crossSectionHandler)
  {
    crossSectionHandler->Clear();
    delete crossSectionHandler;
  }

  if (shellCrossSectionHandler)
    {
      shellCrossSectionHandler->Clear();
      delete shellCrossSectionHandler;
    }

  
  // Reading of data files - all materials are read  
  crossSectionHandler = new G4CrossSectionHandler;
  crossSectionHandler->Clear();
  G4String crossSectionFile = "phot/pe-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  shellCrossSectionHandler = new G4CrossSectionHandler();
  shellCrossSectionHandler->Clear();
  G4String shellCrossSectionFile = "phot/pe-ss-cs-";
  shellCrossSectionHandler->LoadShellData(shellCrossSectionFile);

  fParticleChange = GetParticleChangeForGamma();

  fAtomDeexcitation = G4LossTableManager::Instance()->AtomDeexcitation();
  if(fAtomDeexcitation) { 
    fDeexcitationActive = fAtomDeexcitation->IsFluoActive(); 
  }
  if (verboseLevel > 1) {
    G4cout << "Livermore Polarized PhotoElectric model is initialized " 
	   << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / keV << " keV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePolarizedPhotoElectricModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  G4double cs = crossSectionHandler->FindValue(G4lrint(Z), GammaEnergy);

  if (verboseLevel > 3) {
    G4cout << "G4LivermorePolarizedPhotoElectricModel::ComputeCrossSectionPerAtom()" 
	   << G4endl;
    G4cout << "  E(keV)= " << GammaEnergy/keV << "  Z= " << Z 
	   << "  CrossSection(barn)= "
	   << cs/barn << G4endl; 
  }
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePolarizedPhotoElectricModel::SampleSecondaries(
                                std::vector<G4DynamicParticle*>* fvect,
				const G4MaterialCutsCouple* couple,
				const G4DynamicParticle* aDynamicGamma,
				G4double,
				G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling SampleSecondaries() of G4LivermorePolarizedPhotoElectricModel" 
	   << G4endl;
  }

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();
  
  // kill incident photon  
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  
  // low-energy gamma is absorpted by this process

  if (photonEnergy <= lowEnergyLimit)
    {
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
      return;
    }

  // Select randomly one element in the current material
  //G4cout << "Select random atom Egamma(keV)= " << photonEnergy/keV << G4endl;
  const G4Element* elm = 
    SelectRandomAtom(couple->GetMaterial(),theGamma,photonEnergy);
  G4int Z = G4lrint(elm->GetZ());

  // Select the ionised shell in the current atom according to shell cross sections
  //G4cout << "Select random shell Z= " << Z << G4endl;
  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);
  //G4cout << "Shell index= " << shellIndex 
  //	 << " nShells= " << elm->GetNbOfAtomicShells() << G4endl;
  G4double bindingEnergy;
  const G4AtomicShell* shell = 0;
  if(fDeexcitationActive) {
    G4AtomicShellEnumerator as = G4AtomicShellEnumerator(shellIndex);
    shell = fAtomDeexcitation->GetAtomicShell(Z, as);
    bindingEnergy = shell->BindingEnergy();
  } else {
    G4int nshells = elm->GetNbOfAtomicShells() - 1;
    if(G4int(shellIndex) > nshells) { shellIndex = std::max(0, nshells); }
    bindingEnergy = elm->GetAtomicShell(shellIndex);
  }

  // There may be cases where the binding energy of the selected 
  // shell is > photon energy
  // In such cases do not generate secondaries
  if(photonEnergy < bindingEnergy) {
    fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
    return;
  }
  //G4cout << "Z= " << Z <<  "  shellIndex= " << shellIndex 
  //	 << " Ebind(keV)= " <<  bindingEnergy/keV << G4endl; 


  // Primary outcoming electron
  G4double eKineticEnergy = photonEnergy - bindingEnergy;
  G4double edep = bindingEnergy;

  // Calculate direction of the photoelectron
  G4ThreeVector electronDirection = 
    GetAngularDistribution()->SampleDirection(aDynamicGamma, 
					      eKineticEnergy,
					      shellIndex, 
					      couple->GetMaterial());

  // The electron is created ...
  G4DynamicParticle* electron = new G4DynamicParticle (theElectron,
						       electronDirection,
						       eKineticEnergy);
  fvect->push_back(electron);

  // Sample deexcitation
  if(fDeexcitationActive) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      size_t nbefore = fvect->size();

      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
	for (size_t j=nbefore; j<nafter; ++j) {
	  edep -= ((*fvect)[j])->GetKineticEnergy();
	} 
      }
    }
  }

  // energy balance - excitation energy left
  if(edep > 0.0) {
    fParticleChange->ProposeLocalEnergyDeposit(edep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
