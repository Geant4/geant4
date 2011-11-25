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
// $Id: G4LivermorePhotoElectricModel.cc,v 1.13 2010-12-27 17:45:12 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyPhotoElectric developed by A.Forti and M.G.Pia
//
// History:
// --------
// 15 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - remove GetMeanFreePath method and table
//                  - simplify sampling of deexcitation by using cut in energy
//                  - added protection against numerical problem in energy sampling 
//                  - use G4ElementSelector
// 23 Oct 2009   L Pandola
//                  - atomic deexcitation managed via G4VEmModel::DeexcitationFlag() is 
//                    set as "true" (default would be false)
// 15 Mar 2010   L Pandola
//                  - removed methods to set explicitely fluorescence cuts.
//                  Main cuts from G4ProductionCutsTable are always used
// 26 Dec 2010   V Ivanchenko Load data tables only once to avoid memory leak
// 30 May 2011   A Mantero & V Ivanchenko Migration to model design for deexcitation
// 

#include "G4LivermorePhotoElectricModel.hh"
#include "G4LossTableManager.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4CrossSectionHandler.hh"
#include "G4ShellData.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4VPhotoElectricAngularDistribution.hh"
#include "G4PhotoElectricAngularGeneratorSimple.hh"
#include "G4PhotoElectricAngularGeneratorSauterGavrila.hh"
#include "G4PhotoElectricAngularGeneratorPolarized.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::G4LivermorePhotoElectricModel(
           const G4ParticleDefinition*, const G4String& nam)
  : G4VEmModel(nam),isInitialised(false),crossSectionHandler(0),
    shellCrossSectionHandler(0),fAtomDeexcitation(0),fElectronAngularGenerator(0)
{
  lowEnergyLimit  = 250 * eV; 
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
  fTransitionManager = G4AtomicTransitionManager::Instance();

  // default generator
  fElectronAngularGenerator = 
    new G4PhotoElectricAngularGeneratorSauterGavrila("GEANTSauterGavrilaGenerator");        

  if(verboseLevel>0) {
    G4cout << "Livermore PhotoElectric is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / eV << " eV - "
	   << highEnergyLimit / GeV << " GeV"
	   << G4endl;
  }

  //Mark this model as "applicable" for atomic deexcitation
  SetDeexcitationFlag(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::~G4LivermorePhotoElectricModel()
{  
  delete crossSectionHandler;
  delete shellCrossSectionHandler;
  delete fElectronAngularGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePhotoElectricModel::Initialise(const G4ParticleDefinition*,
					  const G4DataVector&)
{
  if (verboseLevel > 2) {
    G4cout << "Calling G4LivermorePhotoElectricModel::Initialise()" << G4endl;
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

  // Read data tables for all materials  
  crossSectionHandler = new G4CrossSectionHandler();
  G4String crossSectionFile = "phot/pe-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  shellCrossSectionHandler = new G4CrossSectionHandler();
  G4String shellCrossSectionFile = "phot/pe-ss-cs-";
  shellCrossSectionHandler->LoadShellData(shellCrossSectionFile);
  
  //  
  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for Livermore PhotoElectric model" << G4endl;
  }


  if(isInitialised) { return; }
  isInitialised = true;

  fParticleChange = GetParticleChangeForGamma();

  fAtomDeexcitation  = G4LossTableManager::Instance()->AtomDeexcitation();

  if (verboseLevel > 0) { 
    G4cout << "Livermore PhotoElectric model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePhotoElectricModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3) {
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4LivermorePhotoElectricModel" 
	   << G4endl;
  }
  if (GammaEnergy < lowEnergyLimit || GammaEnergy > highEnergyLimit) {
    return 0;
  }
  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePhotoElectricModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
						 const G4MaterialCutsCouple* couple,
						 const G4DynamicParticle* aDynamicGamma,
						 G4double,
						 G4double)
{

  // Fluorescence generated according to:
  // J. Stepanek ,"A program to determine the radiation spectra due to a single atomic
  // subshell ionisation by a particle or due to deexcitation or decay of radionuclides", 
  // Comp. Phys. Comm. 1206 pp 1-1-9 (1997)
 
  if (verboseLevel > 3) {
    G4cout << "Calling SampleSecondaries() of G4LivermorePhotoElectricModel" << G4endl;
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
 
  // Returns the normalized direction of the momentum
  G4ThreeVector photonDirection = aDynamicGamma->GetMomentumDirection(); 

  // Select randomly one element in the current material
  const G4Element* elm = SelectRandomAtom(couple->GetMaterial(),theGamma,photonEnergy);
  G4int Z = (G4int)elm->GetZ();

  // Select the ionised shell in the current atom according to shell cross sections
  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);
  //G4double bindingEnergy = elm->GetAtomicShell(shellIndex);

  // Retrieve the corresponding identifier and binding energy of the selected shell
  const G4AtomicShell* shell = fTransitionManager->Shell(Z,shellIndex);
  G4double bindingEnergy = shell->BindingEnergy();
  G4int shellId = shell->ShellId();

  // Primary outcoming electron
  G4double eKineticEnergy = photonEnergy - bindingEnergy;

  // There may be cases where the binding energy of the selected shell is > photon energy
  // In such cases do not generate secondaries
  if (eKineticEnergy > 0.)
    {
      // Calculate direction of the photoelectron
      G4ThreeVector gammaPolarization = aDynamicGamma->GetPolarization();
      G4ThreeVector electronDirection = 
	fElectronAngularGenerator->GetPhotoElectronDirection(photonDirection,
							     eKineticEnergy,
							     gammaPolarization,
							     shellId);

      // The electron is created ...
      G4DynamicParticle* electron = new G4DynamicParticle (theElectron,
							   electronDirection,
							   eKineticEnergy);
      fvect->push_back(electron);
    }
  else
    {
      bindingEnergy = photonEnergy;
    }


  // Sample deexcitation
  if(fAtomDeexcitation) {
    G4int index = couple->GetIndex();
    if(fAtomDeexcitation->CheckDeexcitationActiveRegion(index)) {
      size_t nbefore = fvect->size();
      fAtomDeexcitation->GenerateParticles(fvect, shell, Z, index);
      size_t nafter = fvect->size();
      if(nafter > nbefore) {
	for (size_t j=nbefore; j<nafter; ++j) {
	  bindingEnergy -= ((*fvect)[j])->GetKineticEnergy();
	} 
      }
    }
  }
  // energy balance
  if(bindingEnergy < 0.0) { bindingEnergy = 0.0; }
  
  // excitation energy left
  fParticleChange->ProposeLocalEnergyDeposit(bindingEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePhotoElectricModel::SetAngularGenerator(G4VPhotoElectricAngularDistribution* dist)
{
  delete fElectronAngularGenerator;
  fElectronAngularGenerator = dist;
  fElectronAngularGenerator->PrintGeneratorInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SetAngularGenerator(const G4String& name)
{
  if (name == "default") 
    {
      delete fElectronAngularGenerator;
      fElectronAngularGenerator = 
	new G4PhotoElectricAngularGeneratorSimple("GEANT4LowEnergySimpleGenerator");
      generatorName = name;
    }
  else if (name == "standard")
    {
      delete fElectronAngularGenerator;
      fElectronAngularGenerator = 
	new G4PhotoElectricAngularGeneratorSauterGavrila("GEANT4SauterGavrilaGenerator");
      generatorName = name;
    }
  else if (name == "polarized")
    {
      delete fElectronAngularGenerator;
      fElectronAngularGenerator = 
	new G4PhotoElectricAngularGeneratorPolarized("GEANT4LowEnergyPolarizedGenerator");
      generatorName = name;
    }
  else
    {
      G4Exception("G4LivermorePhotoElectricModel::SetAngularGenerator",
		    "em1008",FatalException,"generator does not exist");
    }

  fElectronAngularGenerator->PrintGeneratorInformation();
}

