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
// $Id: G4PenelopePhotoElectricModel.cc,v 1.13 2010-11-26 11:51:11 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
//
// History:
// --------
// 08 Oct 2008   L Pandola  Migration from process to model 
// 08 Jan 2009   L Pandola  Check shell index to avoid mismatch between 
//                          the Penelope cross section database and the 
//                          G4AtomicTransitionManager database. It suppresses 
//                          a warning from G4AtomicTransitionManager only. 
//                          Results are unchanged.
// 25 Mar 2009   L Pandola  Small fix to avoid wrong energy-violation warnings
// 17 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - do not apply production threshold on secondaries
// 19 May 2009   L Pandola    Explicitely set to zero pointers deleted in 
//                            Initialise(), since they might be checked later on
// 21 Oct 2009   L Pandola    Remove un-necessary fUseAtomicDeexcitation flag - now managed by
//                            G4VEmModel::DeexcitationFlag()
// 15 Mar 2010   L Pandola    Explicitely initialize Auger to false
//

#include "G4PenelopePhotoElectricModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsTable.hh"
#include "G4ElementTable.hh"
#include "G4Element.hh"
#include "G4CrossSectionHandler.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4VEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4PenelopePhotoElectricModel::G4PenelopePhotoElectricModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),crossSectionHandler(0),
   shellCrossSectionHandler(0)
{
  fIntrinsicLowEnergyLimit = 100.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  verboseLevel= 0;
  // Verbosity scale:
  // 0 = nothing 
  // 1 = warning for energy non-conservation 
  // 2 = details of energy budget
  // 3 = calculation of cross sections, file openings, sampling of atoms
  // 4 = entering in methods

  //by default the model will inkove the atomic deexcitation
  SetDeexcitationFlag(true);  
  ActivateAuger(false);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PenelopePhotoElectricModel::~G4PenelopePhotoElectricModel()
{  
  if (crossSectionHandler) delete crossSectionHandler;
  if (shellCrossSectionHandler) delete shellCrossSectionHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopePhotoElectricModel::Initialise(const G4ParticleDefinition*,
                                       const G4DataVector& )
{
  if (verboseLevel > 3)
    G4cout << "Calling  G4PenelopePhotoElectricModel::Initialise()" << G4endl;
  if (crossSectionHandler)
    {
      crossSectionHandler->Clear();
      delete crossSectionHandler;
      crossSectionHandler = 0;
    }
  if (shellCrossSectionHandler)
    {
      shellCrossSectionHandler->Clear();
      delete shellCrossSectionHandler;
      shellCrossSectionHandler =0;
    }

  //Re-initialize cross section handlers
  crossSectionHandler = new G4CrossSectionHandler();
  crossSectionHandler->Clear();
  G4String crossSectionFile = "penelope/ph-cs-pen-";
  crossSectionHandler->LoadData(crossSectionFile);
  shellCrossSectionHandler = new G4CrossSectionHandler();
  shellCrossSectionHandler->Clear();
  crossSectionFile = "penelope/ph-ss-cs-pen-";
  shellCrossSectionHandler->LoadShellData(crossSectionFile);
  //This is used to retrieve cross section values later on
  G4VEMDataSet* emdata = 
    crossSectionHandler->BuildMeanFreePathForMaterials();
  //The method BuildMeanFreePathForMaterials() is required here only to force 
  //the building of an internal table: the output pointer can be deleted
  delete emdata;

  if (verboseLevel > 2) 
    G4cout << "Loaded cross section files for PenelopePhotoElectric" << G4endl;

  if (verboseLevel > 0) { 
    G4cout << "Penelope Photo-Electric model is initialized " << G4endl
	   << "Energy range: "
	   << LowEnergyLimit() / MeV << " MeV - "
	   << HighEnergyLimit() / GeV << " GeV"
	   << G4endl;
  }

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForGamma();
  isInitialised = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopePhotoElectricModel::ComputeCrossSectionPerAtom(
                                       const G4ParticleDefinition*,
                                             G4double energy,
                                             G4double Z, G4double,
                                             G4double, G4double)
{
  //
  // Penelope model. Use data-driven approach for cross section estimate (and 
  // also shell sampling from a given atom). Data are from the Livermore database
  //  D.E. Cullen et al., Report UCRL-50400 (1989)
  //

  if (verboseLevel > 3)
    G4cout << "Calling ComputeCrossSectionPerAtom() of G4PenelopePhotoElectricModel" << G4endl;

  G4int iZ = (G4int) Z;
  //  if (!crossSectionHandler) // VI: should not be 
  //  {
  //    G4cout << "G4PenelopePhotoElectricModel::ComputeCrossSectionPerAtom" << G4endl;
  //    G4cout << "The cross section handler is not correctly initialized" << G4endl;
  //    G4Exception();
  //  }
  G4double cs = crossSectionHandler->FindValue(iZ,energy);
 
  if (verboseLevel > 2)
    G4cout << "Photoelectric cross section at " << energy/MeV << " MeV for Z=" << Z <<
      " = " << cs/barn << " barn" << G4endl;
  return cs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopePhotoElectricModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
					      const G4MaterialCutsCouple* couple,
					      const G4DynamicParticle* aDynamicGamma,
					      G4double,
					      G4double)
{
  //
  // Photoelectric effect, Penelope model
  //
  // The target atom and the target shell are sampled according to the Livermore 
  // database 
  //  D.E. Cullen et al., Report UCRL-50400 (1989)
  // The angular distribution of the electron in the final state is sampled 
  // according to the Sauter distribution from 
  //  F. Sauter, Ann. Phys. 11 (1931) 454
  // The energy of the final electron is given by the initial photon energy minus 
  // the binding energy. Fluorescence de-excitation is subsequently produced 
  // (to fill the vacancy) according to the general Geant4 G4DeexcitationManager:
  //  J. Stepanek, Comp. Phys. Comm. 1206 pp 1-1-9 (1997)

  if (verboseLevel > 3)
    G4cout << "Calling SamplingSecondaries() of G4PenelopePhotoElectricModel" << G4endl;

  G4double photonEnergy = aDynamicGamma->GetKineticEnergy();

  // always kill primary
  fParticleChange->ProposeTrackStatus(fStopAndKill);
  fParticleChange->SetProposedKineticEnergy(0.);

  if (photonEnergy <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->ProposeLocalEnergyDeposit(photonEnergy);
      return ;
    }

  G4ParticleMomentum photonDirection = aDynamicGamma->GetMomentumDirection();

  // Select randomly one element in the current material
  if (verboseLevel > 2)
    G4cout << "Going to select element in " << couple->GetMaterial()->GetName() << G4endl;
  //use crossSectionHandler instead of G4EmElementSelector because in this case 
  //the dimension of the table is equal to the dimension of the database 
  //(less interpolation errors)
  G4int Z = crossSectionHandler->SelectRandomAtom(couple,photonEnergy);
  if (verboseLevel > 2)
    G4cout << "Selected Z = " << Z << G4endl;

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

  G4double localEnergyDeposit = 0.0;

  // Primary outcoming electron
  G4double eKineticEnergy = photonEnergy - bindingEnergy;
  
  // There may be cases where the binding energy of the selected shell is > photon energy
  // In such cases do not generate secondaries
  if (eKineticEnergy > 0.)
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
      fvect->push_back(electron);
    } 
  else
    {
      bindingEnergy = photonEnergy;
    }
  G4double energyInFluorescence = 0; //testing purposes

  //Now, take care of fluorescence, if required
  if(DeexcitationFlag() && Z > 5) 
    {
      const G4ProductionCutsTable* theCoupleTable=
	G4ProductionCutsTable::GetProductionCutsTable();
      size_t indx = couple->GetIndex();
      G4double cutG = (*(theCoupleTable->GetEnergyCutsVector(0)))[indx];
      G4double cutE = (*(theCoupleTable->GetEnergyCutsVector(1)))[indx];
      
      // Protection to avoid generating photons in the unphysical case of 
      // shell binding energy > photon energy
      if (bindingEnergy > cutG || bindingEnergy > cutE)
	{
	  deexcitationManager.SetCutForSecondaryPhotons(cutG);
	  deexcitationManager.SetCutForAugerElectrons(cutE);
	  std::vector<G4DynamicParticle*>* photonVector = 
	    deexcitationManager.GenerateParticles(Z,shellId); 
	  //Check for secondaries
          if(photonVector) 
	    {
	      for (size_t k=0; k< photonVector->size(); k++)
		{
		  G4DynamicParticle* aPhoton = (*photonVector)[k];
		  if (aPhoton)
		    {
		      G4double itsEnergy = aPhoton->GetKineticEnergy();
		      if (itsEnergy <= bindingEnergy)
			{
			  if(aPhoton->GetDefinition() == G4Gamma::Gamma())
			    energyInFluorescence += itsEnergy;
			  bindingEnergy -= itsEnergy;
			  fvect->push_back(aPhoton);
			}
		      else
			{
			  delete aPhoton;
			  (*photonVector)[k] = 0;
			}
		    }
		}
	      delete photonVector;
	    }
	}
    }
  //Residual energy is deposited locally
  localEnergyDeposit += bindingEnergy;
      
  if (localEnergyDeposit < 0)
    {
      G4cout << "WARNING - "
	     << "G4PenelopePhotoElectric::PostStepDoIt - Negative energy deposit"
	     << G4endl;
      localEnergyDeposit = 0;
    }

  fParticleChange->ProposeLocalEnergyDeposit(localEnergyDeposit);

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4PenelopePhotoElectric" << G4endl;
      G4cout << "Incoming photon energy: " << photonEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      if (eKineticEnergy)
	G4cout << "Outgoing electron " << eKineticEnergy/keV << " keV" << G4endl;
      G4cout << "Fluorescence: " << energyInFluorescence/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << localEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (eKineticEnergy+energyInFluorescence+localEnergyDeposit)/keV << 
	" keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      G4double energyDiff = 
	std::fabs(eKineticEnergy+energyInFluorescence+localEnergyDeposit-photonEnergy);
      if (energyDiff > 0.05*keV)
	G4cout << "Warning from G4PenelopePhotoElectric: problem with energy conservation: " << 
	  (eKineticEnergy+energyInFluorescence+localEnergyDeposit)/keV 
	       << " keV (final) vs. " << 
	  photonEnergy/keV << " keV (initial)" << G4endl;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4PenelopePhotoElectricModel::ActivateAuger(G4bool augerbool)
{
  if (!DeexcitationFlag() && augerbool)
    {
      G4cout << "WARNING - G4PenelopePhotoElectricModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(augerbool);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << augerbool << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopePhotoElectricModel::SampleElectronDirection(G4double energy)
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

