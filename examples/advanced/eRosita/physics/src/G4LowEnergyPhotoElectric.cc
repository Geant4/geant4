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
//
// Author: A. Forti
//         Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// --------
// October 1998 - low energy modifications by Alessandra Forti
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Added EnergySampling method A. Forti
// Modified PostStepDoIt to insert sampling with EPDL97 data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
//   10.04.2000 VL
// - Correcting Fluorescence transition probabilities in order to take into account 
//   non-radiative transitions. No Auger electron simulated yet: energy is locally deposited.
// 17.02.2000 Veronique Lefebure
// - bugs corrected in fluorescence simulation: 
//   . when final use of binding energy: no photon was ever created
//   . no Fluorescence was simulated when the photo-electron energy
//     was below production threshold.
//
// 07-09-99,  if no e- emitted: edep=photon energy, mma
// 24.04.01   V.Ivanchenko     remove RogueWave 
// 12.08.2001 MGP              Revised according to a design iteration
// 16.09.2001 E. Guardincerri  Added fluorescence generation
// 06.10.2001 MGP              Added protection to avoid negative electron energies
//                             when binding energy of selected shell > photon energy
// 18.04.2001 V.Ivanchenko     Fix problem with low energy gammas from fluorescence
//                             MeanFreePath is calculated by crosSectionHandler directly
// 31.05.2002 V.Ivanchenko     Add path of Fluo + Auger cuts to AtomicDeexcitation
// 14.06.2002 V.Ivanchenko     By default do not cheak range of e-
// 21.01.2003 V.Ivanchenko     Cut per region
// 10.05.2004 P.Rodrigues      Changes to accommodate new angular generators
// 20.01.2006 A.Trindade       Changes to accommodate polarized angular generator
//
// --------------------------------------------------------------

#include "G4LowEnergyPhotoElectric.hh"

#include "G4RDVPhotoElectricAngularDistribution.hh"
#include "G4RDPhotoElectricAngularGeneratorSimple.hh"
#include "G4RDPhotoElectricAngularGeneratorSauterGavrila.hh"
#include "G4RDPhotoElectricAngularGeneratorPolarized.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ThreeVector.hh"
#include "G4RDVCrossSectionHandler.hh"
#include "G4RDCrossSectionHandler.hh"
#include "G4RDVEMDataSet.hh"
#include "G4RDCompositeEMDataSet.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLogLogInterpolation.hh"
#include "G4RDVRangeTest.hh"
#include "G4RDRangeNoTest.hh"
#include "G4RDAtomicTransitionManager.hh"
#include "G4RDAtomicShell.hh"
#include "G4ProductionCutsTable.hh"

G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric(const G4String& processName)
  : G4VDiscreteProcess(processName), lowEnergyLimit(250*eV), highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV),
    cutForLowEnergySecondaryPhotons(250.*eV),
    cutForLowEnergySecondaryElectrons(250.*eV)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric()",
                  "OutOfRange", FatalException,
                  "Energy limit outside intrinsic process validity range!");
    }

  crossSectionHandler = new G4RDCrossSectionHandler();
  shellCrossSectionHandler = new G4RDCrossSectionHandler();
  meanFreePathTable = 0;
  rangeTest = new G4RDRangeNoTest;
  generatorName = "geant4.6.2";
  ElectronAngularGenerator = new G4RDPhotoElectricAngularGeneratorSimple("GEANTSimpleGenerator");              // default generator


  if (verboseLevel > 0)
    {
      G4cout << GetProcessName() << " is created " << G4endl
	     << "Energy range: "
	     << lowEnergyLimit / keV << " keV - "
	     << highEnergyLimit / GeV << " GeV"
	     << G4endl;
    }
}

G4LowEnergyPhotoElectric::~G4LowEnergyPhotoElectric()
{
  delete crossSectionHandler;
  delete shellCrossSectionHandler;
  delete meanFreePathTable;
  delete rangeTest;
  delete ElectronAngularGenerator;
}

void G4LowEnergyPhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& )
{

  crossSectionHandler->Clear();
  G4String crossSectionFile = "phot/pe-cs-";
  crossSectionHandler->LoadData(crossSectionFile);

  shellCrossSectionHandler->Clear();
  G4String shellCrossSectionFile = "phot/pe-ss-cs-";
  shellCrossSectionHandler->LoadShellData(shellCrossSectionFile);

  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4LowEnergyPhotoElectric::PostStepDoIt(const G4Track& aTrack,
							  const G4Step& aStep)
{
  // Fluorescence generated according to:
  // J. Stepanek ,"A program to determine the radiation spectra due to a single atomic
  // subshell ionisation by a particle or due to deexcitation or decay of radionuclides", 
  // Comp. Phys. Comm. 1206 pp 1-1-9 (1997)
 
  aParticleChange.Initialize(aTrack);
  
  const G4DynamicParticle* incidentPhoton = aTrack.GetDynamicParticle();
  G4double photonEnergy = incidentPhoton->GetKineticEnergy();
  if (photonEnergy <= lowEnergyLimit)
    {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      aParticleChange.ProposeEnergy(0.);
      aParticleChange.ProposeLocalEnergyDeposit(photonEnergy);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }
 
  G4ThreeVector photonDirection = incidentPhoton->GetMomentumDirection(); // Returns the normalized direction of the momentum

  // Select randomly one element in the current material
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy);

  // Select the ionised shell in the current atom according to shell cross sections
  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);

  // Retrieve the corresponding identifier and binding energy of the selected shell
  const G4RDAtomicTransitionManager* transitionManager = G4RDAtomicTransitionManager::Instance();
  const G4RDAtomicShell* shell = transitionManager->Shell(Z,shellIndex);
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
      // Generate the electron only if with large enough range w.r.t. cuts and safety
      G4double safety = aStep.GetPostStepPoint()->GetSafety();

      if (rangeTest->Escape(G4Electron::Electron(),couple,eKineticEnergy,safety))
	{

	  // Calculate direction of the photoelectron
	  G4ThreeVector gammaPolarization = incidentPhoton->GetPolarization();
	  G4ThreeVector electronDirection = ElectronAngularGenerator->GetPhotoElectronDirection(photonDirection,eKineticEnergy,gammaPolarization,shellId);

	  // The electron is created ...
	  G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(),
							       electronDirection,
							       eKineticEnergy);
	  electronVector.push_back(electron);
	}
      else
	{
	  energyDeposit += eKineticEnergy;
	}
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

  G4int nSecondaries  = nElectrons + nPhotons;
  aParticleChange.SetNumberOfSecondaries(nSecondaries);

  for (G4int l = 0; l<nElectrons; l++ )
    {
      aPhoton = electronVector[l];
      if(aPhoton) {
        aParticleChange.AddSecondary(aPhoton);
      }
    }
  for ( size_t ll = 0; ll < nTotPhotons; ll++)
    {
      aPhoton = (*photonVector)[ll];
      if(aPhoton) {
        aParticleChange.AddSecondary(aPhoton);
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

  // Kill the incident photon
  aParticleChange.ProposeMomentumDirection( 0., 0., 0. );
  aParticleChange.ProposeEnergy( 0. );

  aParticleChange.ProposeLocalEnergyDeposit(energyDeposit);
  aParticleChange.ProposeTrackStatus( fStopAndKill );

  // Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

G4bool G4LowEnergyPhotoElectric::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() );
}

G4double G4LowEnergyPhotoElectric::GetMeanFreePath(const G4Track& track,
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

void G4LowEnergyPhotoElectric::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForLowEnergySecondaryPhotons = cut;
  deexcitationManager.SetCutForSecondaryPhotons(cut);
}

void G4LowEnergyPhotoElectric::SetCutForLowEnSecElectrons(G4double cut)
{
  cutForLowEnergySecondaryElectrons = cut;
  deexcitationManager.SetCutForAugerElectrons(cut);
}

void G4LowEnergyPhotoElectric::ActivateAuger(G4bool val)
{
  deexcitationManager.ActivateAugerElectronProduction(val);
}

void G4LowEnergyPhotoElectric::SetAngularGenerator(G4RDVPhotoElectricAngularDistribution* distribution)
{
  ElectronAngularGenerator = distribution;
  ElectronAngularGenerator->PrintGeneratorInformation();
}

void G4LowEnergyPhotoElectric::SetAngularGenerator(const G4String& name)
{
  if (name == "default") 
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = new G4RDPhotoElectricAngularGeneratorSimple("GEANT4LowEnergySimpleGenerator");
      generatorName = name;
    }
  else if (name == "standard")
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = new G4RDPhotoElectricAngularGeneratorSauterGavrila("GEANT4SauterGavrilaGenerator");
      generatorName = name;
    }
  else if (name == "polarized")
    {
      delete ElectronAngularGenerator;
      ElectronAngularGenerator = new G4RDPhotoElectricAngularGeneratorPolarized("GEANT4LowEnergyPolarizedGenerator");
      generatorName = name;
    }
  else
    {
      G4Exception("G4LowEnergyPhotoElectric::SetAngularGenerator()",
                  "InvalidSetup", FatalException,
                  "Generator does not exist!");
    }

  ElectronAngularGenerator->PrintGeneratorInformation();
}
