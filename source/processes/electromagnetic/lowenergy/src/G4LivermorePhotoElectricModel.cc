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
// Author: Sebastien Inserti
//         30 October 2008
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
// 

#include "G4LivermorePhotoElectricModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::G4LivermorePhotoElectricModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false),meanFreePathTable(0),
 crossSectionHandler(0),shellCrossSectionHandler(0),ElectronAngularGenerator(0)
{
  lowEnergyLimit = 250 * eV; 
  highEnergyLimit = 100 * GeV;
  //  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);
  
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  //Set atomic deexcitation by default
  SetDeexcitationFlag(true);
  ActivateAuger(false);

  // default generator
  ElectronAngularGenerator = 
    new G4PhotoElectricAngularGeneratorSauterGavrila("GEANTSauterGavrilaGenerator");        

  if(verboseLevel>0) {
    G4cout << "Livermore PhotoElectric is constructed " << G4endl
	   << "Energy range: "
	   << lowEnergyLimit / eV << " eV - "
	   << highEnergyLimit / GeV << " GeV"
	   << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::~G4LivermorePhotoElectricModel()
{  
  delete crossSectionHandler;
  delete shellCrossSectionHandler;
  delete ElectronAngularGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePhotoElectricModel::Initialise(const G4ParticleDefinition*,
					  const G4DataVector&)
{
  if (verboseLevel > 3) {
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
  //crossSectionHandler->Clear();
  G4String crossSectionFile = "phot/pe-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  shellCrossSectionHandler = new G4CrossSectionHandler();
  //shellCrossSectionHandler->Clear();
  G4String shellCrossSectionFile = "phot/pe-ss-cs-";
  shellCrossSectionHandler->LoadShellData(shellCrossSectionFile);
  
  //  
  if (verboseLevel > 2) {
    G4cout << "Loaded cross section files for Livermore PhotoElectric model" << G4endl;
  }
  //  InitialiseElementSelectors(particle,cuts);

  if (verboseLevel > 0) { 
    G4cout << "Livermore PhotoElectric model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / eV << " eV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
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
 
  if (verboseLevel > 3)
    G4cout << "Calling SampleSecondaries() of G4LivermorePhotoElectricModel" << G4endl;

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
  //  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy);
  const G4ParticleDefinition* particle =  aDynamicGamma->GetDefinition();
  const G4Element* elm = SelectRandomAtom(couple->GetMaterial(),particle,photonEnergy);
  G4int Z = (G4int)elm->GetZ();

  // Select the ionised shell in the current atom according to shell cross sections
  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);

  // Retrieve the corresponding identifier and binding energy of the selected shell
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  const G4AtomicShell* shell = transitionManager->Shell(Z,shellIndex);
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
	ElectronAngularGenerator->GetPhotoElectronDirection(photonDirection,
							    eKineticEnergy,
							    gammaPolarization,
							    shellId);

      // The electron is created ...
      G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							   electronDirection,
							   eKineticEnergy);
      fvect->push_back(electron);
    }
  else
    {
      bindingEnergy = photonEnergy;
    }

  // deexcitation
  if(DeexcitationFlag() && Z > 5) {

    const G4ProductionCutsTable* theCoupleTable=
      G4ProductionCutsTable::GetProductionCutsTable();

    size_t index = couple->GetIndex();
    G4double cutg = (*(theCoupleTable->GetEnergyCutsVector(0)))[index];
    G4double cute = (*(theCoupleTable->GetEnergyCutsVector(1)))[index];

    // Generation of fluorescence
    // Data in EADL are available only for Z > 5
    // Protection to avoid generating photons in the unphysical case of
    // shell binding energy > photon energy
    if (bindingEnergy > cutg || bindingEnergy > cute)
      {
	G4DynamicParticle* aPhoton;
	deexcitationManager.SetCutForSecondaryPhotons(cutg);
	deexcitationManager.SetCutForAugerElectrons(cute);
 
	std::vector<G4DynamicParticle*>* photonVector = 
	  deexcitationManager.GenerateParticles(Z,shellId);
	size_t nTotPhotons = photonVector->size();
	for (size_t k=0; k<nTotPhotons; k++)
	  {
	    aPhoton = (*photonVector)[k];
	    if (aPhoton)
	      {
		G4double itsEnergy = aPhoton->GetKineticEnergy();
		if (itsEnergy <= bindingEnergy)
		  {
		    // Local energy deposit is given as the sum of the
		    // energies of incident photons minus the energies
		    // of the outcoming fluorescence photons
		    bindingEnergy -= itsEnergy;
		    fvect->push_back(aPhoton);
		  }
		else
		  {
		    // abnormal case of energy non-conservation
		    delete aPhoton;
		    (*photonVector)[k] = 0;
		  }
	      }
	  }
	delete photonVector;
      }
  }
  // excitation energy left
  fParticleChange->ProposeLocalEnergyDeposit(bindingEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::ActivateAuger(G4bool augerbool)
{
  if (!DeexcitationFlag() && augerbool)
    {
      G4cout << "WARNING - G4LivermorePhotoElectricModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(augerbool);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << augerbool << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4LivermorePhotoElectricModel::SetAngularGenerator(G4VPhotoElectricAngularDistribution* dist)
{
  delete ElectronAngularGenerator;
  ElectronAngularGenerator = dist;
  ElectronAngularGenerator->PrintGeneratorInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SetAngularGenerator(const G4String& name)
{
  if (name == "default") 
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = 
	new G4PhotoElectricAngularGeneratorSimple("GEANT4LowEnergySimpleGenerator");
      generatorName = name;
    }
  else if (name == "standard")
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = 
	new G4PhotoElectricAngularGeneratorSauterGavrila("GEANT4SauterGavrilaGenerator");
      generatorName = name;
    }
  else if (name == "polarized")
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = 
	new G4PhotoElectricAngularGeneratorPolarized("GEANT4LowEnergyPolarizedGenerator");
      generatorName = name;
    }
  else
    {
      G4Exception("G4LowEnergyPhotoElectric::SetAngularGenerator - generator does not exist");
    }

  ElectronAngularGenerator->PrintGeneratorInformation();
}

