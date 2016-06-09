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
// $Id: G4LivermoreIonisationModel.cc,v 1.14 2010-12-03 16:03:35 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Luciano Pandola
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
// 
//

#include "G4LivermoreIonisationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProductionCutsTable.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicDeexcitation.hh"
#include "G4AtomicShell.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4CrossSectionHandler.hh"
#include "G4ProcessManager.hh"
#include "G4VEMDataSet.hh"
#include "G4EMDataSet.hh"
#include "G4eIonisationCrossSectionHandler.hh"
#include "G4eIonisationSpectrum.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4ShellVacancy.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4CompositeEMDataSet.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


G4LivermoreIonisationModel::G4LivermoreIonisationModel(const G4ParticleDefinition*,
                                             const G4String& nam)
  :G4VEmModel(nam),isInitialised(false),crossSectionHandler(0),
   energySpectrum(0),shellVacancy(0)
{
  fIntrinsicLowEnergyLimit = 10.0*eV;
  fIntrinsicHighEnergyLimit = 100.0*GeV;
  fNBinEnergyLoss = 360;
  //  SetLowEnergyLimit(fIntrinsicLowEnergyLimit);
  SetHighEnergyLimit(fIntrinsicHighEnergyLimit);
  //
  verboseLevel = 0;
  //By default: use deexcitation, not auger
  SetDeexcitationFlag(true);
  ActivateAuger(false);
  //
  //
  // Notice: the fluorescence along step is generated only if it is 
  // set by the PROCESS (e.g. G4eIonisation) via the command
  // process->ActivateDeexcitation(true);
  //  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LivermoreIonisationModel::~G4LivermoreIonisationModel()
{
  if (energySpectrum) delete energySpectrum;
  if (crossSectionHandler) delete crossSectionHandler;
  if (shellVacancy) delete shellVacancy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreIonisationModel::Initialise(const G4ParticleDefinition* particle,
					    const G4DataVector& cuts)
{
  //Check that the Livermore Ionisation is NOT attached to e+
  if (particle != G4Electron::Electron())
    {
      G4cout << "ERROR: Livermore Ionisation Model is applicable only to electrons" << G4endl;
      G4cout << "It cannot be registered to " << particle->GetParticleName() << G4endl;
      G4Exception();
    }

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

  G4VDataSetAlgorithm* interpolation = new G4SemiLogInterpolation();
  crossSectionHandler = new G4eIonisationCrossSectionHandler(energySpectrum,interpolation,
							     LowEnergyLimit(),HighEnergyLimit(),
							     fNBinEnergyLoss);
  crossSectionHandler->Clear();
  crossSectionHandler->LoadShellData("ioni/ion-ss-cs-");
  //This is used to retrieve cross section values later on
  G4VEMDataSet* emdata = 
    crossSectionHandler->BuildMeanFreePathForMaterials(&cuts);
  //The method BuildMeanFreePathForMaterials() is required here only to force 
  //the building of an internal table: the output pointer can be deleted
  delete emdata;  
 
  //Fluorescence data
  transitionManager = G4AtomicTransitionManager::Instance();
  if (shellVacancy) delete shellVacancy;
  shellVacancy = new G4ShellVacancy();
  InitialiseFluorescence();

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

  if(isInitialised) return;
  fParticleChange = GetParticleChangeForLoss();
  isInitialised = true; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreIonisationModel::MinEnergyCut(const G4ParticleDefinition*,
						  const G4MaterialCutsCouple*)
{
  return 250.*eV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreIonisationModel::ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
								G4double energy,
								G4double Z, G4double,
								G4double cutEnergy, 
								G4double)
{
  G4int iZ = (G4int) Z;
  if (!crossSectionHandler)
    {
      G4cout << "G4LivermoreIonisationModel::ComputeCrossSectionPerAtom" << G4endl;
      G4cout << "The cross section handler is not correctly initialized" << G4endl;
      G4Exception();
      return 0;
    }
  
  //The cut is already included in the crossSectionHandler
  G4double cs = 
    crossSectionHandler->GetCrossSectionAboveThresholdForElement(energy,cutEnergy,iZ);

  if (verboseLevel > 1)
    {
      G4cout << "G4LivermoreIonisationModel " << G4endl;
      G4cout << "Cross section for delta emission > " << cutEnergy/keV << " keV at " <<
	energy/keV << " keV and Z = " << iZ << " --> " << cs/barn << " barn" << G4endl;
    }
  return cs;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4LivermoreIonisationModel::ComputeDEDXPerVolume(const G4Material* material,
							  const G4ParticleDefinition* ,
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
      G4cout << "Stopping power < " << cutEnergy/keV << " keV at " << 
	kineticEnergy/keV << " keV = " << sPower/(keV/mm) << " keV/mm" << G4endl;
    }
  
  return sPower;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreIonisationModel::SampleSecondaries(std::vector<G4DynamicParticle*>* fvect,
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
      return ;
    }

   // Select atom and shell
  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);
  G4int shell = crossSectionHandler->SelectRandomShell(Z, kineticEnergy);
  const G4AtomicShell* atomicShell =
                (G4AtomicTransitionManager::Instance())->Shell(Z, shell);
  G4double bindingEnergy = atomicShell->BindingEnergy();
  G4int shellId = atomicShell->ShellId();

  // Sample delta energy using energy interval for delta-electrons 
  G4double energyMax = 
    std::min(maxE,energySpectrum->MaxEnergyOfSecondaries(kineticEnergy));
  G4double energyDelta = energySpectrum->SampleEnergy(Z, cutE, energyMax,
						      kineticEnergy, shell);

  if (energyDelta == 0.) //nothing happens
    return;

  // Transform to shell potential
  G4double deltaKinE = energyDelta + 2.0*bindingEnergy;
  G4double primaryKinE = kineticEnergy + 2.0*bindingEnergy;

  // sampling of scattering angle neglecting atomic motion
  G4double deltaMom = std::sqrt(deltaKinE*(deltaKinE + 2.0*electron_mass_c2));
  G4double primaryMom = std::sqrt(primaryKinE*(primaryKinE + 2.0*electron_mass_c2));

  G4double cost = deltaKinE * (primaryKinE + 2.0*electron_mass_c2)
                            / (deltaMom * primaryMom);
  if (cost > 1.) cost = 1.;
  G4double sint = std::sqrt(1. - cost*cost);
  G4double phi  = twopi * G4UniformRand();
  G4double dirx = sint * std::cos(phi);
  G4double diry = sint * std::sin(phi);
  G4double dirz = cost;
  
  // Rotate to incident electron direction
  G4ThreeVector primaryDirection = aDynamicParticle->GetMomentumDirection();
  G4ThreeVector deltaDir(dirx,diry,dirz);
  deltaDir.rotateUz(primaryDirection);
  //Updated components
  dirx = deltaDir.x();
  diry = deltaDir.y();
  dirz = deltaDir.z();

  // Take into account atomic motion del is relative momentum of the motion
  // kinetic energy of the motion == bindingEnergy in V.Ivanchenko model
  cost = 2.0*G4UniformRand() - 1.0;
  sint = std::sqrt(1. - cost*cost);
  phi  = twopi * G4UniformRand();
  G4double del = std::sqrt(bindingEnergy *(bindingEnergy + 2.0*electron_mass_c2))
               / deltaMom;
  dirx += del* sint * std::cos(phi);
  diry += del* sint * std::sin(phi);
  dirz += del* cost;

  // Find out new primary electron direction
  G4double finalPx = primaryMom*primaryDirection.x() - deltaMom*dirx;
  G4double finalPy = primaryMom*primaryDirection.y() - deltaMom*diry;
  G4double finalPz = primaryMom*primaryDirection.z() - deltaMom*dirz;

  //Ok, ready to create the delta ray
  G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
  theDeltaRay->SetKineticEnergy(energyDelta);
  G4double norm = 1.0/std::sqrt(dirx*dirx + diry*diry + dirz*dirz);
  dirx *= norm;
  diry *= norm;
  dirz *= norm;
  theDeltaRay->SetMomentumDirection(dirx, diry, dirz);
  theDeltaRay->SetDefinition(G4Electron::Electron());
  fvect->push_back(theDeltaRay);

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
      G4double norm = 1.0/std::sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
      finalPx *= norm;
      finalPy *= norm;
      finalPz *= norm;
      fParticleChange->ProposeMomentumDirection(finalPx, finalPy, finalPz);
    }
  fParticleChange->SetProposedKineticEnergy(finalKinEnergy);

  // deexcitation may be active per G4Region
  if(DeexcitationFlag() && Z > 5) 
    {
      G4ProductionCutsTable* theCoupleTable = 
	G4ProductionCutsTable::GetProductionCutsTable();
      // Retrieve cuts for gammas
      G4double cutG = (*theCoupleTable->GetEnergyCutsVector(0))[couple->GetIndex()];
      if(theEnergyDeposit > cutG || theEnergyDeposit > cutE) 
	{
	  deexcitationManager.SetCutForSecondaryPhotons(cutG);
	  deexcitationManager.SetCutForAugerElectrons(cutE);
	  std::vector<G4DynamicParticle*>* secondaryVector = 
	    deexcitationManager.GenerateParticles(Z, shellId);
	  G4DynamicParticle* aSecondary;

	  if (secondaryVector) 
	    {
	      for (size_t i = 0; i<secondaryVector->size(); i++) 
		{
		  aSecondary = (*secondaryVector)[i];
		  //Check if it is a valid secondary
		  if (aSecondary) 
		    {
		      G4double e = aSecondary->GetKineticEnergy();
		      if (e < theEnergyDeposit) 
			{		      
			  theEnergyDeposit -= e;
			  fvect->push_back(aSecondary);
			  aSecondary = 0;
			  (*secondaryVector)[i]=0;
			} 
		      else 
			{
			  delete aSecondary;
			  (*secondaryVector)[i] = 0;
			}
		    }
		}
	      //secondaryVector = 0; 
	      delete secondaryVector;
	    }
	}
    }

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
      G4cout << "Total final state: " << (finalKinEnergy+energyDelta+bindingEnergy+
					  theEnergyDeposit)/keV << " keV" << G4endl;
      G4cout << "-----------------------------------------------------------" << G4endl;
    }
  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void G4LivermoreIonisationModel::SampleDeexcitationAlongStep(const G4Material* theMaterial,
							     const G4Track& theTrack,
							     G4double& eloss)
{
  //No call if there is no deexcitation along step
  if (!DeexcitationFlag()) return;

  //This method gets the energy loss calculated "Along the Step" and 
  //(including fluctuations) and produces explicit fluorescence/Auger 
  //secondaries. The eloss value is updated.
  G4double energyLossBefore = eloss;

  if (verboseLevel > 2)
    {
      G4cout << "-----------------------------------------------------------" << G4endl;
      G4cout << " SampleDeexcitationAlongStep() from G4LivermoreIonisation" << G4endl;
      G4cout << "Energy loss along step before deexcitation : " << energyLossBefore/keV << 
	" keV" << G4endl;
    }
  G4double incidentEnergy = theTrack.GetDynamicParticle()->GetKineticEnergy();

  G4ProductionCutsTable* theCoupleTable = 
    G4ProductionCutsTable::GetProductionCutsTable();
  const G4MaterialCutsCouple* couple = theTrack.GetMaterialCutsCouple();
  size_t index = couple->GetIndex();
  G4double cutg = (*(theCoupleTable->GetEnergyCutsVector(0)))[index];
  G4double cute = (*(theCoupleTable->GetEnergyCutsVector(1)))[index];
  

  std::vector<G4DynamicParticle*>* deexcitationProducts =
    new std::vector<G4DynamicParticle*>;

  if(eloss > cute || eloss > cutg)
    {
      const G4AtomicTransitionManager* transitionManager =
	G4AtomicTransitionManager::Instance();
      deexcitationManager.SetCutForSecondaryPhotons(cutg);
      deexcitationManager.SetCutForAugerElectrons(cute);
      
      size_t nElements = theMaterial->GetNumberOfElements();
      const G4ElementVector* theElementVector = theMaterial->GetElementVector();

      std::vector<G4DynamicParticle*>* secVector = 0;
      G4DynamicParticle* aSecondary = 0;
      //G4ParticleDefinition* type = 0;
      G4double e;
      G4ThreeVector position;
      G4int shell, shellId;

      // sample secondaries  
      G4double eTot = 0.0;
      std::vector<G4int> n =
	shellVacancy->GenerateNumberOfIonisations(couple,
						  incidentEnergy,eloss);
      for (size_t i=0; i<nElements; i++) 
	{
	  G4int Z = (G4int)((*theElementVector)[i]->GetZ());
	  size_t nVacancies = n[i];
	  G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();
	  if (nVacancies > 0 && Z > 5 && (maxE > cute || maxE > cutg))
	    {     
	      for (size_t j=0; j<nVacancies; j++) 
		{
		  shell = crossSectionHandler->SelectRandomShell(Z, incidentEnergy);
		  shellId = transitionManager->Shell(Z, shell)->ShellId();
		  G4double maxEShell =
		    transitionManager->Shell(Z, shell)->BindingEnergy();
		  if (maxEShell > cute || maxEShell > cutg ) 
		    {
		      secVector = deexcitationManager.GenerateParticles(Z, shellId);
		      if (secVector) 
			{
			  for (size_t l = 0; l<secVector->size(); l++) {
			    aSecondary = (*secVector)[l];
			    if (aSecondary) 
			      {
				e = aSecondary->GetKineticEnergy();
				G4double itsCut = cutg;
				if (aSecondary->GetParticleDefinition() == G4Electron::Electron())
				  itsCut = cute;
				if ( eTot + e <= eloss && e > itsCut )
				  {
				    eTot += e;
				    deexcitationProducts->push_back(aSecondary);
				  } 
				else
				  { 
				    delete aSecondary;
				  }
			      }
			  }
			  delete secVector;
			}
		    }
		}
	    }
	}
    }

  G4double energyLossInFluorescence = 0.0;
  size_t nSecondaries = deexcitationProducts->size();
  if (nSecondaries > 0) 
    {
      //You may have already secondaries produced by SampleSubCutSecondaries()
      //at the process G4VEnergyLossProcess
      G4int secondariesBefore = fParticleChange->GetNumberOfSecondaries();
      fParticleChange->SetNumberOfSecondaries(nSecondaries+secondariesBefore);
      const G4StepPoint* preStep = theTrack.GetStep()->GetPreStepPoint();
      const G4StepPoint* postStep = theTrack.GetStep()->GetPostStepPoint();
      G4ThreeVector r = preStep->GetPosition();
      G4ThreeVector deltaR = postStep->GetPosition();
      deltaR -= r;
      G4double t = preStep->GetGlobalTime();
      G4double deltaT = postStep->GetGlobalTime();
      deltaT -= t;
      G4double time, q;
      G4ThreeVector position;

      for (size_t i=0; i<nSecondaries; i++) 
	{
	  G4DynamicParticle* part = (*deexcitationProducts)[i];
	  if (part) 
	    {
	      G4double eSecondary = part->GetKineticEnergy();
	      eloss -= eSecondary;
	      if (eloss > 0.)
		{
		  q = G4UniformRand();
		  time = deltaT*q + t;
		  position  = deltaR*q;
		  position += r;
		  G4Track* newTrack = new G4Track(part, time, position);
	          energyLossInFluorescence += eSecondary;
		  pParticleChange->AddSecondary(newTrack);
		}
	      else
		{
		  eloss += eSecondary;
		  delete part;
		  part = 0;
		}
	    }
	}
    }
  delete deexcitationProducts;
  
  //Check and verbosities. Ensure energy conservation
  if (verboseLevel > 2)
    {
      G4cout << "Energy loss along step after deexcitation : " << eloss/keV <<  
	" keV" << G4endl;
    }
  if (verboseLevel > 1)
    {
      G4cout << "------------------------------------------------------------------" << G4endl;
      G4cout << "Energy in fluorescence: " << energyLossInFluorescence/keV << " keV" << G4endl;
      G4cout << "Residual energy loss: " << eloss/keV << " keV " << G4endl;
      G4cout << "Total final: " << (energyLossInFluorescence+eloss)/keV << " keV" << G4endl;
      G4cout << "Total initial: " << energyLossBefore/keV << " keV" << G4endl;	
      G4cout << "------------------------------------------------------------------" << G4endl;
    }
  if (verboseLevel > 0)
    {
      if (std::fabs(energyLossBefore-energyLossInFluorescence-eloss)>10*eV)
	{
	  G4cout << "Found energy non-conservation at SampleDeexcitationAlongStep() " << G4endl;
	  G4cout << "Energy in fluorescence: " << energyLossInFluorescence/keV << " keV" << G4endl;
	  G4cout << "Residual energy loss: " << eloss/keV << " keV " << G4endl;
	  G4cout << "Total final: " << (energyLossInFluorescence+eloss)/keV << " keV" << G4endl;
	  G4cout << "Total initial: " << energyLossBefore/keV << " keV" << G4endl;	
	}
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreIonisationModel::InitialiseFluorescence()
{
  G4DataVector* ksi = 0;
  G4DataVector* energyVector = 0;
  size_t binForFluo = fNBinEnergyLoss/10;

  //Used to produce a log-spaced energy grid. To be deleted at the end.
  G4PhysicsLogVector* eVector = new G4PhysicsLogVector(LowEnergyLimit(),HighEnergyLimit(),
						       binForFluo);
  const G4ProductionCutsTable* theCoupleTable=
    G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  // Loop on couples
  for (size_t m=0; m<numOfCouples; m++) 
    {
       const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
       const G4Material* material= couple->GetMaterial();

       const G4ElementVector* theElementVector = material->GetElementVector();
       size_t NumberOfElements = material->GetNumberOfElements() ;
       const G4double* theAtomicNumDensityVector =
	 material->GetAtomicNumDensityVector();

       G4VDataSetAlgorithm* interp = new G4LogLogInterpolation();
       G4VEMDataSet* xsis = new G4CompositeEMDataSet(interp, 1., 1.);
       //loop on elements
       G4double energyCut = (*(theCoupleTable->GetEnergyCutsVector(1)))[m];
       for (size_t iel=0; iel<NumberOfElements; iel++ ) 
	 {
	   G4int iZ = (G4int)((*theElementVector)[iel]->GetZ());
	   energyVector = new G4DataVector();
	   ksi    = new G4DataVector();
	   //Loop on energy
	   for (size_t j = 0; j<binForFluo; j++) 
	     {
	       G4double energy = eVector->GetLowEdgeEnergy(j);
	       G4double cross   = 0.;
	       G4double eAverage= 0.;
	       G4int nShells = transitionManager->NumberOfShells(iZ);

	       for (G4int n=0; n<nShells; n++) 
		 {
		   G4double e = energySpectrum->AverageEnergy(iZ, 0.0,energyCut,
							      energy, n);
		   G4double pro = energySpectrum->Probability(iZ, 0.0,energyCut,
							      energy, n);
		   G4double cs= crossSectionHandler->FindValue(iZ, energy, n);
		   eAverage   += e * cs * theAtomicNumDensityVector[iel];
		   cross      += cs * pro * theAtomicNumDensityVector[iel];
		   if(verboseLevel > 1) 
		     {
		       G4cout << "Z= " << iZ
			    << " shell= " << n
			    << " E(keV)= " << energy/keV
			    << " Eav(keV)= " << e/keV
			    << " pro= " << pro
			    << " cs= " << cs
			    << G4endl;
		     }
		 }
	       
	       G4double coeff = 0.0;
	       if(eAverage > 0.) 
		 {
		   coeff = cross/eAverage;
		   eAverage /= cross;
		 }
	       
	       if(verboseLevel > 1) 
		 {
		   G4cout << "Ksi Coefficient for Z= " << iZ
			  << " E(keV)= " << energy/keV
			  << " Eav(keV)= " << eAverage/keV
			  << " coeff= " << coeff
			  << G4endl;
		 }	       
	       energyVector->push_back(energy);
	       ksi->push_back(coeff);
	     }
	   G4VDataSetAlgorithm* interp1 = new G4LogLogInterpolation();
	   G4VEMDataSet* p = new G4EMDataSet(iZ,energyVector,ksi,interp1,1.,1.);
	   xsis->AddComponent(p);
	 }
       if(verboseLevel>3) xsis->PrintData();
       shellVacancy->AddXsiTable(xsis);
    }
  if (eVector) 
    delete eVector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LivermoreIonisationModel::ActivateAuger(G4bool val)
{
  if (!DeexcitationFlag() && val)
    {
      G4cout << "WARNING - G4LivermoreIonisationModel" << G4endl;
      G4cout << "The use of the Atomic Deexcitation Manager is set to false " << G4endl;
      G4cout << "Therefore, Auger electrons will be not generated anyway" << G4endl;
    }
  deexcitationManager.ActivateAugerElectronProduction(val);
  if (verboseLevel > 1)
    G4cout << "Auger production set to " << val << G4endl;
}

