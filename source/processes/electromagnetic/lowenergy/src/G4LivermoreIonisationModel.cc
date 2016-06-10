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
// $Id: G4LivermoreIonisationModel.cc 93810 2015-11-02 11:27:56Z gcosmo $
//
// Author: Luciano Pandola
//         on base of G4LowEnergyIonisation developed by A.Forti and V.Ivanchenko
//
// History:
// --------
// 12 Jan 2009   L Pandola    Migration from process to model 
// 03 Mar 2009   L Pandola    Bug fix (release memory in the destructor)
// 15 Apr 2009   V Ivanchenko Cleanup initialisation and generation of secondaries:
//                  - apply internal high-energy limit only in constructor 
//                  - do not apply low-energy limit (default is 0)
//                  - simplify sampling of deexcitation by using cut in energy
//                  - set activation of Auger "false"
//                  - remove initialisation of element selectors
// 19 May 2009   L Pandola    Explicitely set to zero pointers deleted in 
//                            Initialise(), since they might be checked later on
// 23 Oct 2009   L Pandola
//                  - atomic deexcitation managed via G4VEmModel::DeexcitationFlag() is
//                    set as "true" (default would be false)
// 12 Oct 2010   L Pandola
//                  - add debugging information about energy in 
//                    SampleDeexcitationAlongStep()
//                  - generate fluorescence SampleDeexcitationAlongStep() only above 
//                    the cuts.
// 01 Jun 2011   V Ivanchenko general cleanup - all old deexcitation code removed
//

#include "G4LivermoreIonisationModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4ParticleChangeForLoss.hh"
#include "G4Electron.hh"
#include "G4CrossSectionHandler.hh"
#include "G4VEMDataSet.hh"
#include "G4eIonisationCrossSectionHandler.hh"
#include "G4eIonisationSpectrum.hh"
#include "G4VEnergySpectrum.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4DeltaAngle.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4LivermoreIonisationModel::G4LivermoreIonisationModel(const G4ParticleDefinition*,
						       const G4String& nam) : 
  G4VEmModel(nam), fParticleChange(0), 
  isInitialised(false),
  crossSectionHandler(0), energySpectrum(0)
{
  fIntrinsicLowEnergyLimit = 12.*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;

  verboseLevel = 0;
  SetAngularDistribution(new G4DeltaAngle());

  transitionManager = G4AtomicTransitionManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreIonisationModel::~G4LivermoreIonisationModel()
{
  delete energySpectrum;
  delete crossSectionHandler;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreIonisationModel::Initialise(const G4ParticleDefinition* particle,
					    const G4DataVector& cuts)
{
  //Check that the Livermore Ionisation is NOT attached to e+
  if (particle != G4Electron::Electron())
    {
      G4Exception("G4LivermoreIonisationModel::Initialise",
		  "em0002",FatalException,
		  "Livermore Ionisation Model is applicable only to electrons");
    }

  transitionManager->Initialise();

  //Read energy spectrum
  if (energySpectrum) 
    {
      delete energySpectrum;
      energySpectrum = 0;
    }
  energySpectrum = new G4eIonisationSpectrum();
  if (verboseLevel > 3)
    G4cout << "G4VEnergySpectrum is initialized" << G4endl;

  //Initialize cross section handler
  if (crossSectionHandler) 
    {
      delete crossSectionHandler;
      crossSectionHandler = 0;
    }

  const size_t nbins = 20;
  G4double emin = LowEnergyLimit();
  G4double emax = HighEnergyLimit();
  G4int ndec = G4int(std::log10(emax/emin) + 0.5);
  if(ndec <= 0) { ndec = 1; }

  G4VDataSetAlgorithm* interpolation = new G4SemiLogInterpolation();
  crossSectionHandler = 
    new G4eIonisationCrossSectionHandler(energySpectrum,interpolation,
					 emin,emax,nbins*ndec);
  crossSectionHandler->Clear();
  crossSectionHandler->LoadShellData("ioni/ion-ss-cs-");
  //This is used to retrieve cross section values later on
  G4VEMDataSet* emdata = 
    crossSectionHandler->BuildMeanFreePathForMaterials(&cuts);
  //The method BuildMeanFreePathForMaterials() is required here only to force 
  //the building of an internal table: the output pointer can be deleted
  delete emdata;  

  if (verboseLevel > 0)
    {
      G4cout << "Livermore Ionisation model is initialized " << G4endl
	     << "Energy range: "
	     << LowEnergyLimit() / keV << " keV - "
	     << HighEnergyLimit() / GeV << " GeV"
	     << G4endl;
    }

  if (verboseLevel > 3)
    {
      G4cout << "Cross section data: " << G4endl; 
      crossSectionHandler->PrintData();
      G4cout << "Parameters: " << G4endl;
      energySpectrum->PrintData();
    }

  if(isInitialised) { return; }
  fParticleChange = GetParticleChangeForLoss();
  isInitialised = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4LivermoreIonisationModel::ComputeCrossSectionPerAtom(
                            const G4ParticleDefinition*,
			    G4double energy,
			    G4double Z, G4double,
			    G4double cutEnergy, 
			    G4double)
{
  G4int iZ = (G4int) Z;
  if (!crossSectionHandler)
    {
      G4Exception("G4LivermoreIonisationModel::ComputeCrossSectionPerAtom",
		  "em1007",FatalException,
		  "The cross section handler is not correctly initialized");
      return 0;
    }
  
  //The cut is already included in the crossSectionHandler
  G4double cs = 
    crossSectionHandler->GetCrossSectionAboveThresholdForElement(energy,
								 cutEnergy,
								 iZ);

  if (verboseLevel > 1)
    {
      G4cout << "G4LivermoreIonisationModel " << G4endl;
      G4cout << "Cross section for delta emission > " 
	     << cutEnergy/keV << " keV at " 
	     << energy/keV << " keV and Z = " << iZ << " --> " 
	     << cs/barn << " barn" << G4endl;
    }
  return cs;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4LivermoreIonisationModel::ComputeDEDXPerVolume(const G4Material* material,
						 const G4ParticleDefinition*,
						 G4double kineticEnergy,
						 G4double cutEnergy)
{
  G4double sPower = 0.0;

  const G4ElementVector* theElementVector = material->GetElementVector();
  size_t NumberOfElements = material->GetNumberOfElements() ;
  const G4double* theAtomicNumDensityVector =
                    material->GetAtomicNumDensityVector();

  // loop for elements in the material
  for (size_t iel=0; iel<NumberOfElements; iel++ ) 
    {
      G4int iZ = (G4int)((*theElementVector)[iel]->GetZ());
      G4int nShells = transitionManager->NumberOfShells(iZ);
      for (G4int n=0; n<nShells; n++) 
	{
	  G4double e = energySpectrum->AverageEnergy(iZ, 0.0,cutEnergy,
						     kineticEnergy, n);
	  G4double cs= crossSectionHandler->FindValue(iZ,kineticEnergy, n);
	  sPower   += e * cs * theAtomicNumDensityVector[iel];
	}
      G4double esp = energySpectrum->Excitation(iZ,kineticEnergy);
      sPower   += esp * theAtomicNumDensityVector[iel];
    }

  if (verboseLevel > 2)
    {
      G4cout << "G4LivermoreIonisationModel " << G4endl;
      G4cout << "Stopping power < " << cutEnergy/keV 
	     << " keV at " << kineticEnergy/keV << " keV = " 
	     << sPower/(keV/mm) << " keV/mm" << G4endl;
    }
  
  return sPower;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreIonisationModel::SampleSecondaries(
                                 std::vector<G4DynamicParticle*>* fvect,
				 const G4MaterialCutsCouple* couple,
				 const G4DynamicParticle* aDynamicParticle,
				 G4double cutE,
				 G4double maxE)
{
  
  G4double kineticEnergy = aDynamicParticle->GetKineticEnergy();

  if (kineticEnergy <= fIntrinsicLowEnergyLimit)
    {
      fParticleChange->SetProposedKineticEnergy(0.);
      fParticleChange->ProposeLocalEnergyDeposit(kineticEnergy);
      return;
    }

   // Select atom and shell
  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);
  G4int shellIndex = crossSectionHandler->SelectRandomShell(Z, kineticEnergy);
  const G4AtomicShell* shell = transitionManager->Shell(Z,shellIndex);
  G4double bindingEnergy = shell->BindingEnergy();

  // Sample delta energy using energy interval for delta-electrons 
  G4double energyMax = 
    std::min(maxE,energySpectrum->MaxEnergyOfSecondaries(kineticEnergy));
  G4double energyDelta = energySpectrum->SampleEnergy(Z, cutE, energyMax,
						      kineticEnergy, shellIndex);

  if (energyDelta == 0.) //nothing happens
    { return; }

  const G4ParticleDefinition* electron = G4Electron::Electron();
  G4DynamicParticle* delta = new G4DynamicParticle(electron, 
    GetAngularDistribution()->SampleDirectionForShell(aDynamicParticle, energyDelta,
						      Z, shellIndex,
                                                      couple->GetMaterial()),
                                                      energyDelta);

  fvect->push_back(delta);

  // Change kinematics of primary particle
  G4ThreeVector direction = aDynamicParticle->GetMomentumDirection();
  G4double totalMomentum = std::sqrt(kineticEnergy*(kineticEnergy + 2*electron_mass_c2));

  G4ThreeVector finalP = totalMomentum*direction - delta->GetMomentum();
  finalP               = finalP.unit();

  //This is the amount of energy available for fluorescence
  G4double theEnergyDeposit = bindingEnergy;

  // fill ParticleChange
  // changed energy and momentum of the actual particle
  G4double finalKinEnergy = kineticEnergy - energyDelta - theEnergyDeposit;
  if(finalKinEnergy < 0.0) 
    {
      theEnergyDeposit += finalKinEnergy;
      finalKinEnergy    = 0.0;
    } 
  else 
    {
      fParticleChange->ProposeMomentumDirection(finalP);
    }
  fParticleChange->SetProposedKineticEnergy(finalKinEnergy);

  if (theEnergyDeposit < 0)
    {
      G4cout <<  "G4LivermoreIonisationModel: Negative energy deposit: "
	     << theEnergyDeposit/eV << " eV" << G4endl;
      theEnergyDeposit = 0.0;
    }

  //Assign local energy deposit
  fParticleChange->ProposeLocalEnergyDeposit(theEnergyDeposit);

  if (verboseLevel > 1)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Energy balance from G4LivermoreIonisation" << G4endl;
      G4cout << "Incoming primary energy: " << kineticEnergy/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << "Outgoing primary energy: " << finalKinEnergy/keV << " keV" << G4endl;
      G4cout << "Delta ray " << energyDelta/keV << " keV" << G4endl;
      G4cout << "Fluorescence: " << (bindingEnergy-theEnergyDeposit)/keV << " keV" << G4endl;
      G4cout << "Local energy deposit " << theEnergyDeposit/keV << " keV" << G4endl;
      G4cout << "Total final state: " << (finalKinEnergy+energyDelta+bindingEnergy)
					  << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

