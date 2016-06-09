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
// $Id: G4LivermorePhotoElectricModel.cc,v 1.1 2008/10/30 14:16:35 sincerti Exp $
// GEANT4 tag $Name: geant4-09-02 $
//

#include "G4LivermorePhotoElectricModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::G4LivermorePhotoElectricModel(const G4ParticleDefinition*,
                                             const G4String& nam)
:G4VEmModel(nam),isInitialised(false)
{
  lowEnergyLimit = 250 * eV; // SI - Could be 10 eV ?
  highEnergyLimit = 100 * GeV;
  SetLowEnergyLimit(lowEnergyLimit);
  SetHighEnergyLimit(highEnergyLimit);
  
  G4cout << "Livermore Compton is constructed " << G4endl
         << "Energy range: "
         << lowEnergyLimit / keV << " keV - "
         << highEnergyLimit / GeV << " GeV"
         << G4endl;
 
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermorePhotoElectricModel::~G4LivermorePhotoElectricModel()
{  
  delete meanFreePathTable;
  delete crossSectionHandler;
  delete shellCrossSectionHandler;
  delete ElectronAngularGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::Initialise(const G4ParticleDefinition* particle,
                                       const G4DataVector& cuts)
{
  if (verboseLevel > 3)
    G4cout << "Calling G4LivermorePhotoElectricModel::Initialise()" << G4endl;

  InitialiseElementSelectors(particle,cuts);

  // Energy limits
  
  if (LowEnergyLimit() < lowEnergyLimit)
  {
      G4cout << "G4LivermorePhotoElectricModel: low energy limit increased from " <<
        LowEnergyLimit()/eV << " eV to " << lowEnergyLimit << " eV" << 
	G4endl;
      SetLowEnergyLimit(lowEnergyLimit);
  }
 
  if (HighEnergyLimit() > highEnergyLimit)
  {
      G4cout << "G4LivermorePhotoElectricModel: high energy limit decreased from " <<
        HighEnergyLimit()/GeV << " GeV to " << highEnergyLimit << " GeV" 
	     << G4endl;
      SetHighEnergyLimit(highEnergyLimit);
  }

  // Read data tables for all materials
  
  crossSectionHandler = new G4CrossSectionHandler();
  crossSectionHandler->Clear();
  G4String crossSectionFile = "phot/pe-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  meanFreePathTable = 0;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();

  shellCrossSectionHandler = new G4CrossSectionHandler();
  shellCrossSectionHandler->Clear();
  G4String shellCrossSectionFile = "phot/pe-ss-cs-";
  shellCrossSectionHandler->LoadShellData(shellCrossSectionFile);
  
  // SI - Buggy default ?
  //generatorName = "geant4.6.2";
  //ElectronAngularGenerator = new G4PhotoElectricAngularGeneratorSimple("GEANTSimpleGenerator");              // default generator

  //
  
  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for Livermore PhotoElectric model" << G4endl;

  G4cout << "Livermore PhotoElectric model is initialized " << G4endl
         << "Energy range: "
         << LowEnergyLimit() / keV << " keV - "
         << HighEnergyLimit() / GeV << " GeV"
         << G4endl;

  if(isInitialised) return;

  if(pParticleChange)
    fParticleChange = reinterpret_cast<G4ParticleChangeForGamma*>(pParticleChange);
  else
    fParticleChange = new G4ParticleChangeForGamma();

  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePhotoElectricModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double GammaEnergy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4LivermorePhotoElectricModel" << G4endl;

  G4double cs = crossSectionHandler->FindValue(G4int(Z), GammaEnergy);
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
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
  
  if (photonEnergy <= lowEnergyLimit)
    {
      fParticleChange->ProposeTrackStatus(fStopAndKill);
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
      // SI - IS THE FOLLOWING RETURN NECESSARY ?
      return ;
    }
 
  G4ThreeVector photonDirection = aDynamicGamma->GetMomentumDirection(); // Returns the normalized direction of the momentum

  // Select randomly one element in the current material
  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy);

  // Select the ionised shell in the current atom according to shell cross sections
  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);

  // Retrieve the corresponding identifier and binding energy of the selected shell
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  const G4AtomicShell* shell = transitionManager->Shell(Z,shellIndex);
  G4double bindingEnergy = shell->BindingEnergy();
  G4int shellId = shell->ShellId();

  // Create lists of pointers to DynamicParticles (photons and electrons)
  // (Is the electron vector necessary? To be checked)
  std::vector<G4DynamicParticle*>* photonVector = 0;
  std::vector<G4DynamicParticle*> electronVector;

  G4double energyDeposit = 0.0;

  // Primary outcoming electron
  G4double eKineticEnergy = photonEnergy - bindingEnergy;

  // There may be cases where the binding energy of the selected shell is > photon energy
  // In such cases do not generate secondaries
  if (eKineticEnergy > 0.)
    {
      // SI - Removed safety
      
      // Generate the electron only if with large enough range w.r.t. cuts and safety
      //G4double safety = aStep.GetPostStepPoint()->GetSafety();

      //if (rangeTest->Escape(G4Electron::Electron(),couple,eKineticEnergy,safety))
	{

	  // Calculate direction of the photoelectron
	  G4ThreeVector gammaPolarization = aDynamicGamma->GetPolarization();
	  G4ThreeVector electronDirection = ElectronAngularGenerator->GetPhotoElectronDirection(photonDirection,eKineticEnergy,gammaPolarization,shellId);

	  // The electron is created ...
	  G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							       electronDirection,
							       eKineticEnergy);
	  electronVector.push_back(electron);
	}
      /*else
	{
	  energyDeposit += eKineticEnergy;
	}*/
    }
  else
    {
      bindingEnergy = photonEnergy;
    }

  G4int nElectrons = electronVector.size();
  size_t nTotPhotons = 0;
  G4int nPhotons=0;
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();

  size_t index = couple->GetIndex();
  G4double cutg = (*(theCoupleTable->GetEnergyCutsVector(0)))[index];
  cutg = std::min(cutForLowEnergySecondaryPhotons,cutg);
  
  G4double cute = (*(theCoupleTable->GetEnergyCutsVector(1)))[index];
  cute = std::min(cutForLowEnergySecondaryPhotons,cute);

  G4DynamicParticle* aPhoton;

  // Generation of fluorescence
  // Data in EADL are available only for Z > 5
  // Protection to avoid generating photons in the unphysical case of
  // shell binding energy > photon energy
  if (Z > 5  && (bindingEnergy > cutg || bindingEnergy > cute))
    {
      photonVector = deexcitationManager.GenerateParticles(Z,shellId);
      nTotPhotons = photonVector->size();
      for (size_t k=0; k<nTotPhotons; k++)
	{
	  aPhoton = (*photonVector)[k];
	  if (aPhoton)
	    {
              G4double itsCut = cutg;
              if(aPhoton->GetDefinition() == G4Electron::Electron()) itsCut = cute;
	      G4double itsEnergy = aPhoton->GetKineticEnergy();

	      if (itsEnergy > itsCut && itsEnergy <= bindingEnergy)
		{
		  nPhotons++;
		  // Local energy deposit is given as the sum of the
		  // energies of incident photons minus the energies
		  // of the outcoming fluorescence photons
		  bindingEnergy -= itsEnergy;

		}
	      else
		{
                  delete aPhoton;
                  (*photonVector)[k] = 0;
                }
	    }
	}
    }

  energyDeposit += bindingEnergy;

  // Final state
  
  for (G4int l = 0; l<nElectrons; l++ )
    {
      aPhoton = electronVector[l];
      if(aPhoton) {
        fvect->push_back(aPhoton);
      }
    }
  for ( size_t ll = 0; ll < nTotPhotons; ll++)
    {
      aPhoton = (*photonVector)[ll];
      if(aPhoton) {
        fvect->push_back(aPhoton);
      }
    }

  delete photonVector;

  if (energyDeposit < 0)
    {
      G4cout << "WARNING - "
	     << "G4LowEnergyPhotoElectric::PostStepDoIt - Negative energy deposit"
	     << G4endl;
      energyDeposit = 0;
    }

  // kill incident photon
  fParticleChange->ProposeMomentumDirection( 0., 0., 0. );
  fParticleChange->SetProposedKineticEnergy(0.);
  fParticleChange->ProposeTrackStatus(fStopAndKill);   
  fParticleChange->ProposeLocalEnergyDeposit(energyDeposit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForLowEnergySecondaryPhotons = cut;
  deexcitationManager.SetCutForSecondaryPhotons(cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SetCutForLowEnSecElectrons(G4double cut)
{
  cutForLowEnergySecondaryElectrons = cut;
  deexcitationManager.SetCutForAugerElectrons(cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::ActivateAuger(G4bool val)
{
  deexcitationManager.ActivateAugerElectronProduction(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SetAngularGenerator(G4VPhotoElectricAngularDistribution* distribution)
{
  ElectronAngularGenerator = distribution;
  ElectronAngularGenerator->PrintGeneratorInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermorePhotoElectricModel::SetAngularGenerator(const G4String& name)
{
  if (name == "default") 
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = new G4PhotoElectricAngularGeneratorSimple("GEANT4LowEnergySimpleGenerator");
      generatorName = name;
    }
  else if (name == "standard")
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = new G4PhotoElectricAngularGeneratorSauterGavrila("GEANT4SauterGavrilaGenerator");
      generatorName = name;
    }
  else if (name == "polarized")
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = new G4PhotoElectricAngularGeneratorPolarized("GEANT4LowEnergyPolarizedGenerator");
      generatorName = name;
    }
  else
    {
      G4Exception("G4LowEnergyPhotoElectric::SetAngularGenerator - generator does not exist");
    }

  ElectronAngularGenerator->PrintGeneratorInformation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermorePhotoElectricModel::GetMeanFreePath(const G4Track& track,
						   G4double, // previousStepSize
					       G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  G4Material* material = track.GetMaterial();
  //  size_t materialIndex = material->GetIndex();

  G4double meanFreePath = DBL_MAX;

  //  if (energy > highEnergyLimit) 
  //    meanFreePath = meanFreePathTable->FindValue(highEnergyLimit,materialIndex);
  //  else if (energy < lowEnergyLimit) meanFreePath = DBL_MAX;
  //  else meanFreePath = meanFreePathTable->FindValue(energy,materialIndex);

  G4double cross = shellCrossSectionHandler->ValueForMaterial(material,energy);
  if(cross > 0.0) meanFreePath = 1.0/cross;

  return meanFreePath;
}

