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
// $Id: G4LowEnergyIonisationVI.cc,v 1.5 2001-10-18 14:15:27 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// --------------------------------------------------------------
//
// File name:     G4LowEnergyIonisationVI
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
// 18.10.01 V.Ivanchenko    Add fluorescence AlongStepDoIt
// 18.10.01 MGP             Revision to improve code quality and consistency with design
//
// --------------------------------------------------------------

#include "G4LowEnergyIonisationVI.hh"
#include "G4eIonisationSpectrum.hh"
#include "G4eIonisationCrossSectionHandler.hh"
#include "G4AtomicTransitionManager.hh"
#include "G4AtomicShell.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4SemiLogInterpolation.hh"
#include "G4LogLogInterpolation.hh"
#include "G4EMDataSet.hh"
#include "G4VEMDataSet.hh"
#include "G4CompositeEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4ShellVacancy.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
 

G4LowEnergyIonisationVI::G4LowEnergyIonisationVI(const G4String& nam)
  : G4eLowEnergyLoss(nam), 
  crossSectionHandler(0),
  theMeanFreePath(0),
  energySpectrum(0),
  shellVacancy(0)
{
  cutForPhotons = 0.;
  cutForElectrons = 0.;
  verboseLevel = 0;
}


G4LowEnergyIonisationVI::~G4LowEnergyIonisationVI()
{
  delete crossSectionHandler;
  delete energySpectrum;
  delete theMeanFreePath;
  delete shellVacancy;
  cutForDelta.clear();
}


void G4LowEnergyIonisationVI::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyIonisation::BuildPhysicsTable start"
           << G4endl;
      }

  cutForDelta.clear();

  // Create and fill IonisationParameters once
  if( energySpectrum != 0 ) delete energySpectrum;
  energySpectrum = new G4eIonisationSpectrum();

  if(verboseLevel > 0) {
    G4cout << "G4VEnergySpectrum is initialized"
           << G4endl;
      }

  // Create and fill G4CrossSectionHandler once

  if ( crossSectionHandler != 0 ) delete crossSectionHandler;
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


void G4LowEnergyIonisationVI::BuildLossTable(
                        const G4ParticleDefinition& aParticleType)
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

  if (shellVacancy != 0) delete shellVacancy;
  shellVacancy = new G4ShellVacancy();
  G4DataVector* ksi = 0;
  G4DataVector* energy = 0;
  G4int binForFluo = totBin/10;

  // Do you really need to create a G4PhysicsLogVector object just for the modest
  // service of calculating the bin coordinate? 
  // Waste of memory and loss of performance.
  // - MGP -
  G4PhysicsLogVector* bVector = new G4PhysicsLogVector(lowKineticEnergy,
		                		       highKineticEnergy,
						       binForFluo);
  G4AtomicTransitionManager* transitionManager = G4AtomicTransitionManager::Instance();
  
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
    G4double tCut = ((G4Electron::Electron())->GetCutsInEnergy())[J];
    if(tCut > highKineticEnergy) tCut = highKineticEnergy;
    cutForDelta[J] = tCut;

    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector = 
                    material->GetAtomicNumDensityVector();
    if(verboseLevel > 1) {
      G4cout << "Energy loss for material # " << J
             << " tCut(keV)= " << tCut/keV
             << G4endl;
      }
      
    // now comes the loop for the kinetic energy values
    for (size_t i = 0; i<totBin; i++) {

      G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
      G4double ionloss = 0.;    
      
      // loop for elements in the material
      for (size_t iel=0; iel<NumberOfElements; iel++ ) {

        G4int Z = (G4int)((*theElementVector)[iel]->GetZ());

	G4int nShells = transitionManager->NumberOfShells(Z);

        for (G4int n=0; n<nShells; n++) {

          G4double e = energySpectrum->AverageEnergy(Z, 0.0, tCut, 
                                                             lowEdgeEnergy, n);
          G4double pro = energySpectrum->Probability(Z, 0.0, tCut, 
                                                             lowEdgeEnergy, n);
          G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy, n);
          ionloss   += e * cs * pro * theAtomicNumDensityVector[iel];
          if(verboseLevel > 1) {
            G4cout << "Z= " << Z
                   << " shell= " << n
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

    // fill data for fluorescence

    G4VDataSetAlgorithm* interp = new G4LogLogInterpolation();
    G4VEMDataSet* xsis = new G4CompositeEMDataSet(interp, 1., 1.);
    for (size_t iel=0; iel<NumberOfElements; iel++ ) {

      G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
      energy = new G4DataVector();
      ksi    = new G4DataVector();

      for (size_t j = 0; j<binForFluo; j++) {

        G4double lowEdgeEnergy = bVector->GetLowEdgeEnergy(j);
        G4double cross   = 0.;
        G4double eAverage= 0.;
	G4int nShells = transitionManager->NumberOfShells(Z);

        for (G4int n=0; n<nShells; n++) {

          G4double e = energySpectrum->AverageEnergy(Z, 0.0, tCut, 
                                                             lowEdgeEnergy, n);
          G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy, n);
          eAverage   += e * cs * theAtomicNumDensityVector[iel];
          cross      += cs * theAtomicNumDensityVector[iel];
	}

        if(eAverage > 0.) cross /= eAverage;
        else              cross  = 0.;

        energy->push_back(lowEdgeEnergy);
        ksi->push_back(cross);
      }
      interp = new G4LogLogInterpolation();
      G4VEMDataSet* set = new G4EMDataSet(Z,energy,ksi,interp,1.,1.);
      xsis->AddComponent(set);
    }
    if(verboseLevel) xsis->PrintData();
    shellVacancy->AddXsiTable(xsis);    
  }
  delete bVector;
}


G4VParticleChange* G4LowEnergyIonisationVI::PostStepDoIt(const G4Track& track,
					                 const G4Step&  step)
{
  // Delta electron production mechanism on base of the model
  // J. Stepanek " A program to determine the radiation spectra due 
  // to a single atomic subshell ionisation by a particle or due to 
  // deexcitation or decay of radionuclides", 
  // Comp. Phys. Comm. 1206 pp 1-19 (1997)

  aParticleChange.Initialize(track);

  const G4Material* material = track.GetMaterial();
  G4double kineticEnergy = track.GetKineticEnergy();

  G4int index = material->GetIndex();
  G4double tCut = cutForDelta[index];

  G4double tmax = energySpectrum->MaxEnergyOfSecondaries(kineticEnergy);

  // Select atom and shell

  G4int Z = crossSectionHandler->SelectRandomAtom(material, kineticEnergy);
  G4int shell = crossSectionHandler->SelectRandomShell(Z, kineticEnergy);
  const G4AtomicShell* atomicShell = (G4AtomicTransitionManager::Instance())->Shell(Z, shell);
  G4double bindingEnergy = atomicShell->BindingEnergy();
  G4int shellId = atomicShell->ShellId();

  // Sample delta energy

  G4double tDelta = energySpectrum->SampleEnergy(Z, tCut, tmax, 
                                                 kineticEnergy, shell);

  if(tDelta == 0.0) 
    return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);

  // Transform to shell potential
  G4double deltaKinE = tDelta + 2.0*bindingEnergy; 
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
  G4DynamicParticle* theDeltaRay = new G4DynamicParticle();
  theDeltaRay->SetKineticEnergy(tDelta);
  theDeltaRay->SetMomentumDirection(dirx, diry, dirz); 
  theDeltaRay->SetDefinition(G4Electron::Electron());
     
  G4double theEnergyDeposit = bindingEnergy;


  // Generation of Fluorescence and Auger
  size_t nSecondaries = 0;
  size_t totalNumber  = 1;
  G4std::vector<G4DynamicParticle*>* secondaryVector = 0;
  G4DynamicParticle* aSecondary = 0;
  G4ParticleDefinition* type = 0;
  G4double em = G4std::min(bindingEnergy, theEnergyDeposit);
 
  // Fluorescence data start from element 6

  if (Z > 5 && (em >= cutForPhotons || em >= cutForElectrons)) {

    secondaryVector = deexcitationManager.GenerateParticles(Z, shellId);

    if (secondaryVector != 0) {

      nSecondaries = secondaryVector->size();
      for (size_t i = 0; i<nSecondaries; i++) {

        aSecondary = (*secondaryVector)[i];
        if (aSecondary) {
	  
          G4double e = aSecondary->GetKineticEnergy();
          type = aSecondary->GetDefinition();
          if (e < theEnergyDeposit && 
                ((type == G4Gamma::Gamma() && e > cutForPhotons ) || 
                 (type == G4Electron::Electron() && e > cutForElectrons ))) {

             theEnergyDeposit -= e;
             totalNumber++;

	  } else {

             delete aSecondary;
             (*secondaryVector)[i] = 0;
	  }
	}
      }
    }
  }
      
  // Save delta-electrons

  aParticleChange.SetNumberOfSecondaries(totalNumber);
  aParticleChange.AddSecondary(theDeltaRay);

  // Save Fluorescence and Auger

  if (secondaryVector) {

    for (size_t l = 0; l < nSecondaries; l++) {

      aSecondary = (*secondaryVector)[l];
      if(aSecondary) aParticleChange.AddSecondary(aSecondary); 
    }
    delete secondaryVector;
  }
     
  // fill ParticleChange 
  // changed energy and momentum of the actual particle

  G4double finalKinEnergy = kineticEnergy - tDelta - theEnergyDeposit;
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


void G4LowEnergyIonisationVI::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database,";
  comments += "Gamma energy sampled from a parametrised formula.";
  comments += "Implementation of the continuous dE/dx part.";  
  comments += "\n At present it can be used for electrons ";
  comments += " in the energy range [250eV,100GeV]";
  comments += "\n the process must work with G4LowEnergyBremsstrahlung";
  
  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}         

G4bool G4LowEnergyIonisationVI::IsApplicable(const G4ParticleDefinition& particle)
{
   return ( (&particle == G4Electron::Electron()) );
}

/*
G4std::vector<G4Track*>* G4LowEnergyIonisationVI::SecondariesAlongStep(
                                            const G4Step& step,
                                                  G4double edep) 
{
  if(edep < cutForPhotons && edep < cutForElectrons) return 0;

  G4AtomicTransitionManager* transitionManager = 
                             G4AtomicTransitionManager::Instance();

  // material

  const G4Material* material = (step.GetTrack())->GetMaterial();
  size_t nElements = material->GetNumberOfElements();
  const G4ElementVector* theElementVector = material->GetElementVector();
  G4bool stop = true;

  for (size_t j=0; j<nElements; j++) {

    G4int Z = (G4int)((*theElementVector)[j]->GetZ());
    G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

    if (Z>5 && (maxE>cutForPhotons || maxE>cutForElectrons) ) {
      stop = false;
      break;
    }
  }

  if(stop) return 0;

  // Step parameters

  const G4ParticleDefinition* part = (step.GetTrack())->GetDefinition();
  const G4StepPoint* preStep = step.GetPreStepPoint();
  const G4StepPoint* postStep = step.GetPostStepPoint();
  G4ThreeVector r = preStep->GetPosition();
  G4ThreeVector deltaR = postStep->GetPosition();
  deltaR -= r;
  G4double t = preStep->GetGlobalTime();
  G4double deltaT = postStep->GetGlobalTime();
  deltaT -= t;
  G4double inEnergy = preStep->GetKineticEnergy();

  // create vector of tracks of secondary particles

  G4std::vector<G4Track*>* vectorTracks = new G4std::vector<G4Track*>;
  G4Track* track = 0;
  G4std::vector<G4DynamicParticle*>* secVector = 0;
  G4DynamicParticle* aSecondary = 0;
  G4ParticleDefinition* type = 0;
  G4double e, time, q;
  G4ThreeVector position;
  G4int shell, shellId;

  // sample secondaries
 
  G4std::vector<G4int> n = shellVacancy->GenerateNumberOfIonisations(material,
                                         inEnergy, edep);
  for (size_t i=0; i<nElements; i++) {

    size_t nVacancies = n[i];
    G4int Z = (G4int)((*theElementVector)[i]->GetZ());
    G4double maxE = transitionManager->Shell(Z, 0)->BindingEnergy();

    if (nVacancies && Z>5 && (maxE>cutForPhotons || maxE>cutForElectrons) ) {

      for(size_t j=0; j<nVacancies; j++) {
     
        q = G4UniformRand();
        shell = crossSectionHandler->SelectRandomShell(Z, inEnergy - q*edep);
        shellId = transitionManager->Shell(Z, shell)->ShellId();
        G4double maxE = transitionManager->Shell(Z, shell)->BindingEnergy();
 
        if (maxE>cutForPhotons || maxE>cutForElectrons ) {
          secVector = deexcitationManager.GenerateParticles(Z, shellId);
	} else {
          secVector = 0;
	}

        if (secVector) {	

          time = deltaT*q + t;
          position  = deltaR*q;
          position += r;

          for (size_t l = 0; l<secVector->size(); l++) {

            aSecondary = (*secVector)[l];
            if(aSecondary) {
	  
              e = aSecondary->GetKineticEnergy();
              type = aSecondary->GetDefinition();
              if ( (type == G4Gamma::Gamma() && e > cutForPhotons ) || 
                   (type == G4Electron::Electron() && e > cutForElectrons) ) {
                     
                     track = new G4Track(aSecondary, time, position);
                     vectorTracks->push_back(track);

	      } else {
                     delete aSecondary;
	      }
	    }
	  }
          delete secVector;
	} 
      }
    }
  }
  if(vectorTracks->size()==0) {
    delete vectorTracks;
    return 0;
  }
  return vectorTracks;
}
*/

G4double G4LowEnergyIonisationVI::GetMeanFreePath(const G4Track& track,
						G4double previousStepSize,
						G4ForceCondition* cond)
{
   *cond = NotForced;
   G4int index = (track.GetMaterial())->GetIndex();
   const G4VEMDataSet* data = theMeanFreePath->GetComponent(index);
   G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
   return meanFreePath; 
} 

void G4LowEnergyIonisationVI::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForPhotons = cut;
}   

void G4LowEnergyIonisationVI::SetCutForLowEnSecElectrons(G4double cut)
{
  cutForElectrons = cut;
}

