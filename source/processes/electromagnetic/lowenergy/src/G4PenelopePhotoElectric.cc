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
// $Id: G4PenelopePhotoElectric.cc,v 1.16 2009-06-11 15:47:08 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: L. Pandola
//
// History:
// --------
// January 2003 - Created
// 12 Feb 2003   MG Pia     Migration to "cuts per region"
// 10 Mar 2003 V.Ivanchenko   Remome CutPerMaterial warning
// 31 May 2005  L. Pandola  Added Sauter formula for the sampling of 
//                          the electron direction
// 08 Jan 2009  L. Pandola  Check shell index to avoid mismatch between 
//                          the Penelope cross section database and the 
//                          G4AtomicTransitionManager database. It suppresses 
//                          a warning from G4AtomicTransitionManager only. 
//                          Results are unchanged.
// --------------------------------------------------------------

#include "G4PenelopePhotoElectric.hh"

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
#include "G4RangeNoTest.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"

G4PenelopePhotoElectric::G4PenelopePhotoElectric(const G4String& processName)
  : G4VDiscreteProcess(processName), lowEnergyLimit(250*eV), highEnergyLimit(100*GeV),
    intrinsicLowEnergyLimit(10*eV),
    intrinsicHighEnergyLimit(100*GeV),
    cutForLowEnergySecondaryPhotons(250.*eV),
    cutForLowEnergySecondaryElectrons(250.*eV)
{
  if (lowEnergyLimit < intrinsicLowEnergyLimit || 
      highEnergyLimit > intrinsicHighEnergyLimit)
    {
      G4Exception("G4PenelopePhotoElectric::G4PenelopePhotoElectric - energy limit outside intrinsic process validity range");
    }

  crossSectionHandler = new G4CrossSectionHandler();
  shellCrossSectionHandler = new G4CrossSectionHandler();
  meanFreePathTable = 0;
  rangeTest = new G4RangeNoTest;

  if (verboseLevel > 0) 
    {
      G4cout << GetProcessName() << " is created " << G4endl
	     << "Energy range: " 
	     << lowEnergyLimit / keV << " keV - "
	     << highEnergyLimit / GeV << " GeV" 
	     << G4endl;
    }

   G4cout << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "   The class G4PenelopePhotoElectric is NOT SUPPORTED ANYMORE. " << G4endl;
   G4cout << "   It will be REMOVED with the next major release of Geant4. " << G4endl;
   G4cout << "   Please consult: https://twiki.cern.ch/twiki/bin/view/Geant4/LoweProcesses" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << "*******************************************************************************" << G4endl;
   G4cout << G4endl;

}

G4PenelopePhotoElectric::~G4PenelopePhotoElectric()
{
  delete crossSectionHandler;
  delete shellCrossSectionHandler;
  delete meanFreePathTable;
  delete rangeTest;
}

void G4PenelopePhotoElectric::BuildPhysicsTable(const G4ParticleDefinition& )
{
  
  crossSectionHandler->Clear();
  G4String crossSectionFile = "penelope/ph-cs-pen-";
  crossSectionHandler->LoadData(crossSectionFile);

  shellCrossSectionHandler->Clear();
  G4String shellCrossSectionFile = "penelope/ph-ss-cs-pen-";
  shellCrossSectionHandler->LoadShellData(shellCrossSectionFile);

  delete meanFreePathTable;
  meanFreePathTable = crossSectionHandler->BuildMeanFreePathForMaterials();
}

G4VParticleChange* G4PenelopePhotoElectric::PostStepDoIt(const G4Track& aTrack,
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
 
  G4ParticleMomentum photonDirection = incidentPhoton->GetMomentumDirection();
   
  // Select randomly one element in the current material
  const G4MaterialCutsCouple* couple = aTrack.GetMaterialCutsCouple();
  
  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy);

  // Select the ionised shell in the current atom according to shell cross sections
  size_t shellIndex = shellCrossSectionHandler->SelectRandomShell(Z,photonEnergy);

  // Retrieve the corresponding identifier and binding energy of the selected shell
  const G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();

  //The number of shell cross section possibly reported in the Penelope database 
  //might be different from the number of shells in the G4AtomicTransitionManager
  //(namely, Penelope may contain more shell, especially for very light elements).
  //In order to avoid a warning message from the G4AtomicTransitionManager, I 
  //add this protection. Results are anyway changed, because when G4AtomicTransitionManager
  //has a shellID>maxID, it sets the shellID to the last valid shell. 
  size_t numberOfShells = (size_t) transitionManager->NumberOfShells(Z);
  if (shellIndex >= numberOfShells)
    shellIndex = numberOfShells-1;

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
      // Generate the electron only if with large enough range w.r.t. cuts and safety
      G4double safety = aStep.GetPostStepPoint()->GetSafety();

      if (rangeTest->Escape(G4Electron::Electron(),couple,eKineticEnergy,safety))
	{
	  // The electron is created
	  // Direction sampled from the Sauter distribution
	  G4double cosTheta = SampleElectronDirection(eKineticEnergy);
	  G4double sinTheta = std::sqrt(1-cosTheta*cosTheta);
	  G4double phi = twopi * G4UniformRand() ;
	  G4double dirx = sinTheta * std::cos(phi);
	  G4double diry = sinTheta * std::sin(phi);
	  G4double dirz = cosTheta ;
	  G4ThreeVector electronDirection(dirx,diry,dirz); //electron direction
	  electronDirection.rotateUz(photonDirection);

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
	     << "G4PenelopePhotoElectric::PostStepDoIt - Negative energy deposit"
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

G4bool G4PenelopePhotoElectric::IsApplicable(const G4ParticleDefinition& particle)
{
  return ( &particle == G4Gamma::Gamma() ); 
}

G4double G4PenelopePhotoElectric::GetMeanFreePath(const G4Track& track, 
						  G4double, // previousStepSize
						  G4ForceCondition*)
{
  const G4DynamicParticle* photon = track.GetDynamicParticle();
  G4double energy = photon->GetKineticEnergy();
  G4Material* material = track.GetMaterial();

  G4double meanFreePath = DBL_MAX;

  G4double cross = shellCrossSectionHandler->ValueForMaterial(material,energy);
  if(cross > 0.0) meanFreePath = 1.0/cross;

  return meanFreePath;
}

void G4PenelopePhotoElectric::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForLowEnergySecondaryPhotons = cut;
  deexcitationManager.SetCutForSecondaryPhotons(cut);
}

void G4PenelopePhotoElectric::SetCutForLowEnSecElectrons(G4double cut)
{
  cutForLowEnergySecondaryElectrons = cut;
  deexcitationManager.SetCutForAugerElectrons(cut);
}

void G4PenelopePhotoElectric::ActivateAuger(G4bool val)
{
  deexcitationManager.ActivateAugerElectronProduction(val);
}


G4double G4PenelopePhotoElectric::SampleElectronDirection(G4double energy)
{
  G4double costheta = 1.0;
  if (energy>1*GeV) return costheta;

  //1) initialize energy-dependent variables
  // Variable naming according to Eq. (2.24) of Penelope Manual 
  // (pag. 44)
  G4double gamma = 1.0 + energy/electron_mass_c2;
  G4double gamma2 = gamma*gamma;
  G4double beta = std::sqrt((gamma2-1.0)/gamma2);
  
  // ac corresponds to "A" of Eq. (2.31) 
  // 
  G4double ac = (1.0/beta) - 1.0;
  G4double a1 = 0.5*beta*gamma*(gamma-1.0)*(gamma-2.0);
  G4double a2 = ac + 2.0;
  G4double gtmax = 2.0*(a1 + 1.0/ac);

  G4double tsam = 0;
  G4double gtr = 0;
  
  //2) sampling. Eq. (2.31) of Penelope Manual
  // tsam = 1-std::cos(theta)
  // gtr = rejection function according to Eq. (2.28)
  do{
    G4double rand = G4UniformRand();
    tsam = 2.0*ac * (2.0*rand + a2*std::sqrt(rand)) / (a2*a2 - 4.0*rand);
    gtr = (2.0 - tsam) * (a1 + 1.0/(ac+tsam));
  }while(G4UniformRand()*gtmax > gtr);
  costheta = 1.0-tsam;
  return costheta;
}





