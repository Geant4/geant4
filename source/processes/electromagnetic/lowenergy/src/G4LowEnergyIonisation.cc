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
// $Id: G4LowEnergyIonisation.cc,v 1.64 2001-10-10 17:37:56 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// --------------------------------------------------------------
//
// File name:     G4LowEnergyIonisation
//
// Author:        Alessandra Forti
// 
// Creation date: March 1999
//
// Modifications:
// - 11.04.2000 VL
//   Changing use of float and G4float casts to G4double casts 
//   because of problems with optimisation (bug ?)
//   10.04.2000 VL
// - Correcting Fluorescence transition probabilities in order to take into account 
//   non-radiative transitions. No Auger electron simulated yet: energy is locally deposited.
//   10.04.2000 VL
// - Correction of incident electron final momentum direction
//   07.04.2000 VL+LU
// - First implementation of continuous energy loss
//   22.03.2000 VL
// - 1 bug corrected in SelectRandomAtom method (units)
//   17.02.2000 Veronique Lefebure
// - 5 bugs corrected: 
//   *in Fluorescence, 2 bugs affecting 
//   . localEnergyDeposition and
//   . number of emitted photons that was then always 1 less
//   *in EnergySampling method: 
//   . expon = Parms[13]+1; (instead of uncorrect -1)
//   . rejection /= Parms[6];(instead of uncorrect Parms[7])
//   . Parms[6] is apparently corrupted in the data file (often = 0)  
//     -->Compute normalisation into local variable rejectionMax
//     and use rejectionMax  in stead of Parms[6]
//
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Added EnergySampling method A. Forti
// Modified PostStepDoIt to insert sampling with EEDL data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
// 20.09.00 V.Ivanchenko update fluctuations 
// 24.04.01 V.Ivanchenko remove RogueWave 
// 22.05.01 V.Ivanchenko update calculation of delta-ray kinematic + 
//                       clean up the code 
// 02.08.01 V.Ivanchenko fix energy conservation for small steps 
// 18.08.01 V.Ivanchenko fix energy conservation for pathalogical delta-energy
// 01.10.01 E. Guardincerri Replaced fluorescence generation in PostStepDoIt
//                          according to design iteration
// 04.10.01 MGP             Minor clean-up in the fluo section, removal of
//                          compilation warnings and extra protection to
//                          prevent from accessing a null pointer                                               
// 29.09.01 V.Ivanchenko    revision based on design iteration
// 10.10.01 MGP             Revision to improve code quality and consistency with design
//
// --------------------------------------------------------------

#include "G4LowEnergyIonisation.hh"
#include "G4eIonisationSpectrum.hh"
#include "G4eIonisationCrossSectionHandler.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
 

G4LowEnergyIonisation::G4LowEnergyIonisation(const G4String& nam)
  : G4eLowEnergyLoss(nam), 
  crossSectionHandler(0),
  theMeanFreePath(0),
  energySpectrum(0)
{
  cutForPhotons = 0.;
  cutForElectrons = 0.;
  verboseLevel = 0;
}


G4LowEnergyIonisation::~G4LowEnergyIonisation()
{
  delete crossSectionHandler;
  delete energySpectrum;
  delete theMeanFreePath;
  cutForDelta.clear();
}


void G4LowEnergyIonisation::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyIonisation::BuildPhysicsTable start"
           << G4endl;
      }

  cutForDelta.clear();

  // Create and fill IonisationParameters once
  if( energySpectrum ) delete energySpectrum;
  energySpectrum = new G4eIonisationSpectrum();

  if(verboseLevel > 0) {
    G4cout << "G4VEnergySpectrum is initialized"
           << G4endl;
      }

  // Create and fill G4CrossSectionHandler once

  if( crossSectionHandler ) delete crossSectionHandler;
  G4VDataSetAlgorithm* interpolation = new G4SemiLogInterpolation();
  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  G4int    totBin = GetNbinEloss();
  crossSectionHandler = new G4eIonisationCrossSectionHandler(energySpectrum, 
							     interpolation,
							     lowKineticEnergy, 
							     highKineticEnergy,
							     totBin);
  crossSectionHandler->LoadShellData("ioni/ion-ss-cs-");

  if (verboseLevel > 0) {
    G4cout << GetProcessName() 
           << " is created; Cross section data: " 
           << G4endl;
    crossSectionHandler->PrintData();
    G4cout << "Parameters: " 
           << G4endl;
    energySpectrum->PrintData();
  }

  // Build loss table for IonisationIV

  BuildLossTable(aParticleType);

  if(verboseLevel > 0) {
    G4cout << "The loss table is built"
           << G4endl;
      }

  if (&aParticleType==G4Electron::Electron()) {

    RecorderOfElectronProcess[CounterOfElectronProcess] = (*this).theLossTable;
    CounterOfElectronProcess++;
    PrintInfoDefinition();  

  } else {

    RecorderOfPositronProcess[CounterOfPositronProcess] = (*this).theLossTable;
    CounterOfPositronProcess++;
  }

  // Build mean free path data using cuts values

  if( theMeanFreePath ) delete theMeanFreePath;
  theMeanFreePath = crossSectionHandler->
                    BuildMeanFreePathForMaterials(&cutForDelta);

  if(verboseLevel > 0) {
    G4cout << "The MeanFreePath table is built"
           << G4endl;
      }

  // Build common DEDX table for all ionisation processes
 
  BuildDEDXTable(aParticleType);

  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyIonisation::BuildPhysicsTable end"
           << G4endl;
      }
 
}


void G4LowEnergyIonisation::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  // Build table for energy loss due to soft brems
  // the tables are built for *MATERIALS* binning is taken from LowEnergyLoss

  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  size_t   totBin = GetNbinEloss();
 
  //  create table
  
  if (theLossTable) { 
      theLossTable->clearAndDestroy();
      delete theLossTable;
  }
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const size_t numOfMaterials = G4Material::GetNumberOfMaterials();
  theLossTable = new G4PhysicsTable(numOfMaterials);
  
  // Clean up the vector of cuts

  cutForDelta.resize(numOfMaterials);

  // Loop for materials
  
  for (size_t J=0; J<numOfMaterials; J++) {
    
    // create physics vector and fill it
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(lowKineticEnergy,
		        				 highKineticEnergy,
							 totBin);

    // get material parameters needed for the energy loss calculation
    const G4Material* material= (*theMaterialTable)[J];

    // the cut cannot be below lowest limit
    G4double tcut = ((G4Electron::Electron())->GetCutsInEnergy())[J];
    if(tcut > highKineticEnergy) tcut = highKineticEnergy;
    cutForDelta[J] = tcut;

    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector = 
                    material->GetAtomicNumDensityVector();
    if(verboseLevel > 1) {
      G4cout << "Energy loss for material # " << J
             << " tcut(keV)= " << tcut/keV
             << G4endl;
      }
      
    // now comes the loop for the kinetic energy values
    for (size_t i = 0; i<totBin; i++) {

      G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
      G4double ionloss = 0.;    
      
      // loop for elements in the material
      for (size_t iel=0; iel<NumberOfElements; iel++ ) {

        G4int Z = (G4int)((*theElementVector)[iel]->GetZ());

	// G4int nShells = crossSectionHandler->NumberOfComponents(Z);
	// - MGP - modified: it is not the responsibility of G4VCrossSectionHandler to provide nShells
	G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
	G4int nShells = transitionManager->NumberOfShells(Z);

        for (G4int n=0; n<nShells; n++) {

          G4double e = energySpectrum->AverageEnergy(Z, 0.0, tcut, lowEdgeEnergy, n);
          G4double pro = energySpectrum->Probability(Z, 0.0, tcut, lowEdgeEnergy, n);
          G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy, n);
          ionloss   += e * cs * pro * theAtomicNumDensityVector[iel];
          if(verboseLevel > 2) {
            G4cout << "Z= " << Z
                   << " E(keV)= " << lowEdgeEnergy/keV
                   << " Eav(keV)= " << e/keV
                   << " pro= " << pro
                   << " cs= " << cs
	           << " loss= " << ionloss
                   << G4endl;
          }
        }
      }	
      aVector->PutValue(i,ionloss);
    }
    theLossTable->insert(aVector);
  }
}


inline G4VParticleChange* G4LowEnergyIonisation::PostStepDoIt(const G4Track& track,
							      const G4Step&  step)
{
  // Delta electron production mechanism on base of the model
  // J. Stepanek " A program to determine the radiation spectra due 
  // to a single atomic subshell ionisation by a particle or due to 
  // deexcitation or decay of radionuclides", 
  // Comp. Phys. Comm. 1206 pp 1-19 (1997)

  aParticleChange.Initialize(track);

  const G4Material* mat = track.GetMaterial();
  G4double kineticEnergy = track.GetKineticEnergy();

  G4int index = mat->GetIndex();
  G4double tcut = cutForDelta[index];

  G4double tmax = energySpectrum->MaxEnergyOfSecondaries(kineticEnergy);

  G4int Z = crossSectionHandler->SelectRandomAtom(mat, kineticEnergy);
  G4int shell = crossSectionHandler->SelectRandomShell(Z, kineticEnergy);
  G4double bindingEnergy = (G4AtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();
  G4double tdel = energySpectrum->SampleEnergy(Z, tcut, tmax, kineticEnergy, shell);

  if(tdel == 0.0) 
    return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);

  // Transform to shell potential
  G4double deltaKinE = tdel + 2.0*bindingEnergy; 
  G4double primaryKinE = kineticEnergy + 2.0*bindingEnergy;   

  // sampling of scattering angle neglecting atomic motion
  G4double deltaMom = sqrt(deltaKinE*(deltaKinE + 2.0*electron_mass_c2));
  G4double primaryMom = sqrt(primaryKinE*(primaryKinE + 2.0*electron_mass_c2));
   
  G4double cost = deltaKinE * (primaryKinE + 2.0*electron_mass_c2)
                            / (deltaMom * primaryMom);

  if (cost > 1.) cost = 1.;
  G4double sint = sqrt(1. - cost*cost);
  G4double phi  = twopi * G4UniformRand(); 
  G4double dirx = sint * cos(phi);
  G4double diry = sint * sin(phi);
  G4double dirz = cost;

  // Rotate to incident electron direction
  G4ThreeVector primaryDirection = track.GetMomentumDirection();
  G4ThreeVector deltaDir(dirx,diry,dirz);
  deltaDir.rotateUz(primaryDirection);
  dirx = deltaDir.x();
  diry = deltaDir.y();
  dirz = deltaDir.z();

  // Take into account atomic motion del is relative momentum of the motion
  // kinetic energy of the motion == bindingEnergy in V.Ivanchenko model

  cost = 2.0*G4UniformRand() - 1.0;
  sint = sqrt(1. - cost*cost);
  phi  = twopi * G4UniformRand(); 
  G4double del = sqrt(bindingEnergy *(bindingEnergy + 2.0*electron_mass_c2))
               / deltaMom;
  dirx += del* sint * cos(phi);
  diry += del* sint * sin(phi);
  dirz += del* cost;

  // Find out new primary electron direction
  G4double finalPx = primaryMom*primaryDirection.x() - deltaMom*dirx; 
  G4double finalPy = primaryMom*primaryDirection.y() - deltaMom*diry; 
  G4double finalPz = primaryMom*primaryDirection.z() - deltaMom*dirz; 

  // create G4DynamicParticle object for delta ray
  aParticleChange.SetNumberOfSecondaries(1);
  G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
  theDeltaRay->SetKineticEnergy(tdel);
  theDeltaRay->SetMomentumDirection(dirx, diry, dirz); 
  theDeltaRay->SetDefinition(G4Electron::Electron());
  aParticleChange.AddSecondary(theDeltaRay);
     
  G4double theEnergyDeposit = bindingEnergy;

  // Fluorescence should be implemented here
  // Fluorescence data start from element 6

  if(thePrimShVec.size() != 0) thePrimShVec.clear();
  thePrimShVec.push_back(thePrimaryShell);

  size_t nElectrons = 1;
  size_t nTotPhotons = 0;
  size_t nPhotons = 0;
  
  // Generation of fluorescence
  if (Z > 5)
    {
      photonVector = deexcitationManager.GenerateParticles(Z,shell);
      if (photonVector != 0)
	{
	  nTotPhotons = photonVector->size();
	  for (size_t k=0; k<nTotPhotons; k++)
	    {
	      G4DynamicParticle* aPhoton = (*photonVector)[k];
	      if (aPhoton != 0)
		{
		  G4double itsKineticEnergy = aPhoton->GetKineticEnergy();
		  G4double eDepositTmp = theEnergyDeposit - itsKineticEnergy;
		  if (itsKineticEnergy >= cutForPhotons && eDepositTmp > 0.)
		    {
		      nPhotons++;
		      // Local energy deposit is given as the sum of the 
		      // energies of incident photons minus the energies
		      // of the outcoming fluorescence photons
		      theEnergyDeposit -= itsKineticEnergy;
		    }
		  else
		    {
		      // The current photon would be below threshold,
		      // or it would cause a negative energy deposit
		      delete aPhoton;
		    }
		}
	    }
	}
    }
      
  size_t nSecondaries  = nElectrons + nPhotons;
  
  aParticleChange.SetNumberOfSecondaries(nSecondaries);
      
  for (size_t l = 0; l < nPhotons; l++) 
    {
      aParticleChange.AddSecondary((*photonVector)[l]); 
    }
  
  if (photonVector != 0)
    {
      delete photonVector;
    }

  // fill ParticleChange 
  // changed energy and momentum of the actual particle

  G4double finalKinEnergy = kineticEnergy - tdel - theEnergyDeposit;
  if(finalKinEnergy < 0.0) {
    theEnergyDeposit += finalKinEnergy;
    finalKinEnergy    = 0.0;
  }
  G4double norm = 1.0/sqrt(finalPx*finalPx+finalPy*finalPy+finalPz*finalPz);
  finalPx *= norm;
  finalPy *= norm;
  finalPz *= norm;

  aParticleChange.SetMomentumChange(finalPx, finalPy, finalPz);
  aParticleChange.SetEnergyChange(finalKinEnergy);
  if(theEnergyDeposit < 0.) {
    G4cout << "G4LowEnergyIonisation: Negative energy deposit: " 
           << theEnergyDeposit/eV << " eV" << G4endl;
    theEnergyDeposit = 0.0;
  }
  aParticleChange.SetLocalEnergyDeposit(theEnergyDeposit);

  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


void G4LowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database,";
  comments += "Gamma energy sampled from a parametrised formula.";
  comments += "Implementation of the continuous dE/dx part.";  
  comments += "\n At present it can be used for electrons ";
  comments += " in the energy range [250eV,100GeV]";
  comments += "\n the process must work with G4LowEnergyBremsstrahlung";
  
  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}         

G4bool G4LowEnergyIonisation::IsApplicable(const G4ParticleDefinition& particle)
{
   return ( (&particle == G4Electron::Electron()) );
}


G4double G4LowEnergyIonisation::GetMeanFreePath(const G4Track& track,
						G4double previousStepSize,
						G4ForceCondition* cond)
{
   *cond = NotForced;
   G4int index = (track.GetMaterial())->GetIndex();
   const G4VEMDataSet* data = theMeanFreePath->GetComponent(index);
   G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
   return meanFreePath; 
} 

void G4LowEnergyIonisation::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForPhotons = cut;
}   

void G4LowEnergyIonisation::SetCutForLowEnSecElectrons(G4double cut)
{
  cutForElectrons = cut;
}

