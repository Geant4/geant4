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
//
// $Id: G4LowEnergyPhotoElectric.cc,v 1.43 2002-04-18 13:41:31 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
//                                  
// --------------------------------------------------------------

#include "G4LowEnergyPhotoElectric.hh"

#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ForceCondition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4DynamicParticle.hh"
#include "G4VParticleChange.hh"
#include "G4ThreeVector.hh"
#include "G4VCrossSectionHandler.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VRangeTest.hh"
#include "G4RangeTest.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"

#include "G4CutsPerMaterialWarning.hh"

G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric(const G4String& processName)
  : G4VDiscreteProcess(processName), lowEnergyLimit(250*eV), highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV),
    cutForLowEnergySecondaryPhotons(250.*eV)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4LowEnergyPhotoElectric::G4LowEnergyPhotoElectric - energy limit outside intrinsic process validity range");
    }

  crossSectionHandler = new G4CrossSectionHandler();
  shellCrossSectionHandler = new G4CrossSectionHandler();
  meanFreePathTable = 0;
  rangeTest = new G4RangeTest;

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
}

void G4LowEnergyPhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& photon)
{
  
  G4CutsPerMaterialWarning warning;
  warning.PrintWarning(&photon);

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
      aParticleChange.SetStatusChange(fStopAndKill);
      aParticleChange.SetEnergyChange(0.);
      aParticleChange.SetLocalEnergyDeposit(photonEnergy);
      return G4VDiscreteProcess::PostStepDoIt(aTrack,aStep);
    }
 
  G4ParticleMomentum photonDirection = incidentPhoton->GetMomentumDirection();
   
  // Select randomly one element in the current material
  G4Material* material = aTrack.GetMaterial();
  G4int Z = crossSectionHandler->SelectRandomAtom(material,photonEnergy);

  // Select the ionised shell in the current atom according to shell cross sections
  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);

  // Retrieve the corresponding identifier and binding energy of the selected shell
  G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  const G4AtomicShell* shell = transitionManager->Shell(Z,shellIndex);
  G4double bindingEnergy = shell->BindingEnergy();
  G4int shellId = shell->ShellId();

  // Create lists of pointers to DynamicParticles (photons and electrons)
  // (Is the electron vector necessary? To be checked)
  G4std::vector<G4DynamicParticle*>* photonVector = 0;
  G4std::vector<G4DynamicParticle*> electronVector;

  G4double energyDeposit = 0.0;

  // Primary outcoming electron
  G4double eKineticEnergy = photonEnergy - bindingEnergy;

  // There may be cases where the binding energy of the selected shell is > photon energy
  // In such cases do not generate secondaries
  if (eKineticEnergy > 0.)
    {
      // Generate the electron only if with large enough range w.r.t. cuts and safety
      G4double safety = aStep.GetPostStepPoint()->GetSafety();

      if (rangeTest->Escape(G4Electron::Electron(),material,eKineticEnergy,safety))
	{
	  // The electron is created in the direction of the incident photon ...  
	  G4DynamicParticle* electron = new G4DynamicParticle (G4Electron::Electron(), 
							       photonDirection, 
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
  G4double cut = G4std::min(cutForLowEnergySecondaryPhotons,
                            G4Gamma::Gamma()->GetEnergyThreshold(material));
  G4DynamicParticle* aPhoton;  

  // Generation of fluorescence
  // Data in EADL are available only for Z > 5
  // Protection to avoid generating photons in the unphysical case of 
  // shell binding energy > photon energy
  if (Z > 5  && bindingEnergy > cut)
    {
      photonVector = deexcitationManager.GenerateParticles(Z,shellId); 
      nTotPhotons = photonVector->size();
      for (size_t k=0; k<nTotPhotons; k++)
	{
	  aPhoton = (*photonVector)[k];
	  if (aPhoton)
	    {
	      G4double itsEnergy = aPhoton->GetKineticEnergy();
	      if (itsEnergy > cut && itsEnergy <= bindingEnergy)
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

  for (size_t l = 0; l<nElectrons; l++ )
    {
      aPhoton = electronVector[l];
      if(aPhoton) {
        aParticleChange.AddSecondary(aPhoton);
      }
    }
  for (l = 0; l < nTotPhotons; l++) 
    {
      aPhoton = (*photonVector)[l];
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
  aParticleChange.SetMomentumChange( 0., 0., 0. );
  aParticleChange.SetEnergyChange( 0. );
  
  aParticleChange.SetLocalEnergyDeposit(energyDeposit);  
  aParticleChange.SetStatusChange( fStopAndKill ); 
  
  // Reset NbOfInteractionLengthLeft and return aParticleChange
  return G4VDiscreteProcess::PostStepDoIt( aTrack, aStep );
}

G4bool G4LowEnergyPhotoElectric::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4LowEnergyPhotoElectric::GetMeanFreePath(const G4Track& track, 
					       G4double previousStepSize, 
					       G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  G4Material* material = track.GetMaterial();
  size_t materialIndex = material->GetIndex();

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
}







