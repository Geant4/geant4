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
// --------------------------------------------------------------
//      GEANT 4 class implementation file
//      CERN Geneva Switzerland
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 27 September 2001
//
// Modifications: 
//
// Class Description: 
//
// Bremsstrahlung process based on the model developed  
// by Alessandra Forti, 1999, and Veronique Lefebure, 2000 
//
// -------------------------------------------------------------------
//
// Class description:
// Low Energy electromagnetic process,
// Further documentation available from http://www.ge.infn.it/geant4/lowE
//
// Class Description: End 
//
// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LowEnergyIonisation.hh"
#include "G4eIonisationElectronSpectrum.hh"
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
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LowEnergyIonisation::G4LowEnergyIonisation(const G4String& nam)
  : G4eLowEnergyLoss(nam), 
    theParam(0),
    crossSectionHandler(0),
    theMeanFreePath(0)
{
  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4LowEnergyIonisation::~G4LowEnergyIonisation()
{
  delete crossSectionHandler;
  delete theParam;
  delete theMeanFreePath;
  cutForDelta.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::BuildPhysicsTable(
                         const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyIonisation::BuildPhysicsTable start"
           << G4endl;
      }

  cutForDelta.clear();

  // Create and fill IonisationParameters once
  if( theParam ) delete theParam;
  theParam = new G4eIonisationElectronSpectrum();

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
  crossSectionHandler = new 
           G4eIonisationCrossSectionHandler(theParam, interpolation,
                        lowKineticEnergy, highKineticEnergy, totBin);
  crossSectionHandler->LoadShellData("ioni/ion-ss-cs-");

  if (verboseLevel > 0) {
    G4cout << GetProcessName() 
           << " is created; Cross section data: " 
           << G4endl;
    crossSectionHandler->PrintData();
    G4cout << "Parameters: " 
           << G4endl;
    theParam->PrintData();
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::BuildLossTable(
                          const G4ParticleDefinition& aParticleType)
{
  // Build table for energy loss due to soft brems
  // the tables are built for *MATERIALS* binning is taken from LowEnergyLoss

  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  G4int    totBin = GetNbinEloss();
 
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
    G4int NumberOfElements = material->GetNumberOfElements() ;
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
        G4int nShells = crossSectionHandler->NumberOfComponents(Z);

        for (G4int n=0; n<nShells; n++) {

          G4double e = theParam->AverageEnergy(Z, 0.0, tcut, lowEdgeEnergy, n);
          G4double pro = theParam->Probability(Z, 0.0, tcut, lowEdgeEnergy, n);
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4LowEnergyIonisation::PostStepDoIt(
                                                 const G4Track& track,
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

  G4double tmax = theParam->MaxEnergyOfSecondaries(kineticEnergy);

  G4int Z = crossSectionHandler->SelectRandomAtom(mat, kineticEnergy);
  G4int shell = crossSectionHandler->SelectRandomShell(Z, kineticEnergy);
  G4double bindingEnergy = (G4AtomicTransitionManager::Instance())->
                           Shell(Z, shell)->BindingEnergy();
  G4double tdel = theParam->SampleEnergy(Z, tcut, tmax, kineticEnergy, shell);

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyIonisation::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database,";
           comments += "Gamma energy sampled from a parametrised formula.";
           comments += "Implementation of the continuous dE/dx part.";  
           comments += "\n At present it can be used for electrons ";
           comments += " in the energy range [250eV,100GeV]";
           comments += 
  "\n the process must work with G4LowEnergyBremsstrahlung";
                     
	   G4cout << G4endl << GetProcessName() << ":  " << comments<<G4endl;

}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



