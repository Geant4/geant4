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
// File name:     G4LowEnergyBremsstrahlung
//
// Author:        V.Ivanchenko (Vladimir.Ivantchenko@cern.ch)
// 
// Creation date: 27 September 2001
//
// Modifications: 
//
// --------------------------------------------------------------
//
// Class Description: 
//
// Bremsstrahlung process based on the model developed  
// by Alessandra Forti, 1999, and Veronique Lefebure, 2000 
//
// --------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4BremsstrahlungElectronSpectrum.hh"
#include "G4BremsstrahlungCrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4LowEnergyBremsstrahlung::G4LowEnergyBremsstrahlung(const G4String& nam)
  : G4eLowEnergyLoss(nam), 
    theBR(0),
    crossSectionHandler(0),
    theMeanFreePath(0)
{
  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4LowEnergyBremsstrahlung::~G4LowEnergyBremsstrahlung()
{
  delete crossSectionHandler;
  delete theBR;
  delete theMeanFreePath;
  cutForSecondaryPhotons.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlung::BuildPhysicsTable(
                            const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyBremsstrahlung::BuildPhysicsTable start"
           << G4endl;
      }

  cutForSecondaryPhotons.clear();

  // Create and fill BremsstrahlungParameters once
  if( theBR ) delete theBR;
  theBR = new G4BremsstrahlungElectronSpectrum();

  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyBremsstrahlungSpectrum is initialized"
           << G4endl;
      }

  // Create and fill G4CrossSectionHandler once

  if( crossSectionHandler ) delete crossSectionHandler;
  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();
  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  G4int    totBin = GetNbinEloss();
  crossSectionHandler = new 
           G4BremsstrahlungCrossSectionHandler(theBR, interpolation,
                       lowKineticEnergy, highKineticEnergy, totBin);
  crossSectionHandler->LoadShellData("brem/br-cs-");

  if (verboseLevel > 0) {
    G4cout << GetProcessName() 
           << " is created; Cross section data: " 
           << G4endl;
    crossSectionHandler->PrintData();
    G4cout << "Parameters: " 
           << G4endl;
    theBR->PrintData();
  }

  // Build loss table for Bremsstrahlung

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
                    BuildMeanFreePathForMaterials(&cutForSecondaryPhotons);

  if(verboseLevel > 0) {
    G4cout << "The MeanFreePath table is built"
           << G4endl;
      }

  // Build common DEDX table for all ionisation processes
 
  BuildDEDXTable(aParticleType);

  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyBremsstrahlung::BuildPhysicsTable end"
           << G4endl;
      }
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlung::BuildLossTable(
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

  cutForSecondaryPhotons.resize(numOfMaterials);

  // Loop for materials
  
  for (size_t J=0; J<numOfMaterials; J++) {
    
    // create physics vector and fill it
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(lowKineticEnergy,
		        				 highKineticEnergy,
							 totBin);

    // get material parameters needed for the energy loss calculation
    const G4Material* material= (*theMaterialTable)[J];

    // the cut cannot be below lowest limit
    G4double tcut = G4std::min(highKineticEnergy,
                             ((G4Gamma::Gamma())->GetCutsInEnergy())[J]);
    cutForSecondaryPhotons[J] = tcut;

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
        G4double e = theBR->AverageEnergy(Z, 0.0, tcut, lowEdgeEnergy); 
        G4double pro = theBR->Probability(Z, 0.0, tcut, lowEdgeEnergy); 
        G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy);
        ionloss   += e * cs * pro * theAtomicNumDensityVector[iel];
        if(verboseLevel > 1) {
          G4cout << "Z= " << Z
                 << "; tcut(keV)= " << tcut/keV
                 << "; E(keV)= " << lowEdgeEnergy/keV
                 << "; Eav(keV)= " << e/keV
                 << "; pro= " << pro
                 << "; cs= " << cs
		 << "; loss= " << ionloss
                 << G4endl;
        }
      }	
      aVector->PutValue(i,ionloss);
    }
    theLossTable->insert(aVector);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4VParticleChange* G4LowEnergyBremsstrahlung::PostStepDoIt(
                                                 const G4Track& track,
                                                 const G4Step&  step)
{
  aParticleChange.Initialize(track);

  const G4Material* mat = track.GetMaterial();
  G4double kineticEnergy = track.GetKineticEnergy();
  G4int index = mat->GetIndex();
  G4double tcut = cutForSecondaryPhotons[index];

  // Control limits
  if(tcut >= kineticEnergy) 
     return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);

  G4int Z = crossSectionHandler->SelectRandomAtom(mat, kineticEnergy);

  G4double tgam = theBR->SampleEnergy(Z, tcut, kineticEnergy, kineticEnergy);

  // Sample gamma angle (Z - axis along the parent particle).
  // Universal distribution suggested by L. Urban (Geant3 manual (1993) 
  // Phys211) derived from Tsai distribution (Rev Mod Phys 49,421(1977))

  G4double totalEnergy = kineticEnergy + electron_mass_c2;   

  const G4double a1 = 0.625, a2 = 3.*a1, d = 27.;
  G4double u = - log(G4UniformRand()*G4UniformRand());

  if (9./(9.+d) > G4UniformRand()) u /= a1;
  else                             u /= a2;
    
  G4double theta = u*electron_mass_c2/totalEnergy;
  G4double phi   = twopi * G4UniformRand();
  G4double dirz  = cos(theta);
  G4double sint  = sqrt(1. - dirz*dirz);
  G4double dirx  = sint*cos(phi);
  G4double diry  = sint*sin(phi); 
    
  G4ThreeVector gamDirection (dirx, diry, dirz);
  G4ThreeVector elecDirection = track.GetMomentumDirection();
    
  gamDirection.rotateUz(elecDirection);   
  
  //
  // Update the incident particle 
  //
    
  G4double finalEnergy = kineticEnergy - tgam;  
    
  // Kinematic problem
  if (finalEnergy < 0.) {
    tgam += finalEnergy;
    finalEnergy = 0.0;
  }

  G4double mom = sqrt((totalEnergy + electron_mass_c2)*kineticEnergy);

  G4double finalX = mom*elecDirection.x() - tgam*gamDirection.x();
  G4double finalY = mom*elecDirection.y() - tgam*gamDirection.y();
  G4double finalZ = mom*elecDirection.z() - tgam*gamDirection.z();
      
  aParticleChange.SetNumberOfSecondaries(1);
  aParticleChange.SetMomentumChange(finalX, finalY, finalZ);
  aParticleChange.SetEnergyChange( finalEnergy );

  // create G4DynamicParticle object for the gamma 
  G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
						    gamDirection, tgam);
  aParticleChange.AddSecondary(aGamma); 

  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlung::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database,";
           comments += "Gamma energy sampled from a parametrised formula.";
           comments += "Implementation of the continuous dE/dx part.";  
           comments += "\n At present it can be used for electrons ";
           comments += " in the energy range [250eV,100GeV]";
           comments += 
  "\n the process must work with G4LowEnergyIonisation";
                     
	   G4cout << G4endl << GetProcessName() << ":  " << comments<<G4endl;

}         

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



