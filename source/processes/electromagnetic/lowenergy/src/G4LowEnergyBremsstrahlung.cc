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
//      ------------ G4LowEnergyBremsstrahlung: 
//                   by Alessandra Forti, March 1999
//
// **************************************************************
// 
// 18.04.2000 V.L.
// - First implementation of continuous energy loss.
// 17.02.2000 Veronique Lefebure
//  - correct bug : the gamma energy was not deposited when the gamma was 
//    not produced when its energy was < cutForLowEnergySecondaryPhotons
//
// Added Livermore data table construction methods A. Forti
// Modified BuildMeanFreePath to read new data tables A. Forti
// Modified PostStepDoIt to insert sampling with with EEDL data A. Forti
// Added SelectRandomAtom A. Forti
// Added map of the elements A. Forti
// 20.09.00 update printout V.Ivanchenko
// 24.04.01 V.Ivanchenko remove RogueWave 
// 21.08.01 V.Ivanchenko new design 
// --------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyBremsstrahlungGen.hh"
#include "G4CrossSectionHandler.hh"
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
    shellCrossSectionHandler(0),
    theMeanFreePathData(0),
    theGenerator(0)
{
  tmax.clear();
  cutForLowEnergySecondaryPhotons.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.... 

G4LowEnergyBremsstrahlung::~G4LowEnergyBremsstrahlung()
{
  delete shellCrossSectionHandler;
  delete theGenerator;
  delete theMeanFreePathData;
  tmax.clear();
  cutForLowEnergySecondaryPhotons.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LowEnergyBremsstrahlung::BuildPhysicsTable(
                            const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyBremsstrahlung::BuildPhysicsTable start"
           << G4endl;
      }

  // Create and fill BremsstrahlungParameters once
  if( !theGenerator ) theGenerator = new G4LowEnergyBremsstrahlungGen();
  theGenerator->SetVerbose(verboseLevel);
  theGenerator->Clear();
  theGenerator->Initialize();

  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyBremsstrahlungGen is initialized"
           << G4endl;
      }

  // Create and fill G4CrossSectionHandler once

  if( !shellCrossSectionHandler ) {
    G4VDataSetAlgorithm* shellInterpolation = new G4LogLogInterpolation();
    shellCrossSectionHandler = new G4CrossSectionHandler(shellInterpolation);
    shellCrossSectionHandler->SetVerbose(verboseLevel);
    shellCrossSectionHandler->SetSecondaryGenerator(theGenerator);
    shellCrossSectionHandler->LoadShellData("brem/br-cs-");
  }
  shellCrossSectionHandler->SetVerbose(verboseLevel);

  if (verboseLevel > 0) {
    G4cout << GetProcessName() << " is created " << G4endl;
    shellCrossSectionHandler->PrintData();
    theGenerator->PrintData();
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

  if( theMeanFreePathData ) delete theMeanFreePathData;
  theMeanFreePathData = shellCrossSectionHandler->
     BuildMeanFreePathForMaterials(cutForLowEnergySecondaryPhotons, tmax);

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
  const size_t numOfMaterials = theMaterialTable->length();
  theLossTable = new G4PhysicsTable(numOfMaterials);
  
  // Clean up the vector of cuts

  cutForLowEnergySecondaryPhotons.resize(numOfMaterials);

  // Loop for materials
  
  for (size_t J=0; J<numOfMaterials; J++) {
    
    // create physics vector and fill it
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(lowKineticEnergy,
		        				 highKineticEnergy,
							 totBin);

    // get material parameters needed for the energy loss calculation
    const G4Material* material= (*theMaterialTable)[J];

    // the cut cannot be below lowest limit
    G4double tcut = G4std::max(theGenerator->MinSecondaryEnergy(material),
                              (G4Gamma::Gamma())->GetCutsInEnergy()[J]);
    cutForLowEnergySecondaryPhotons[J] = tcut;

    if(tcut == theGenerator->MinSecondaryEnergy(material) 
            || verboseLevel > 0) {
      G4cout<<"G4LowEnergyBremsstrahlung: gamma Tcut = "
            << tcut/keV
            << " keV for material " << material->GetName()
	    << G4endl;
    }

    const G4ElementVector* theElementVector = material->GetElementVector();
    G4int NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector = 
                    material->GetAtomicNumDensityVector();
      
    // now comes the loop for the kinetic energy values
    for (size_t i = 0; i<totBin; i++) {

      G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
      G4double ionloss = 0.;    
      
      // loop for elements in the material
      for (size_t iel=0; iel<NumberOfElements; iel++ ) {
        G4double Z = (*theElementVector)(iel)->GetZ();
        G4double e = theGenerator->AverageEnergy(Z, lowEdgeEnergy, tcut);
        G4double cs= shellCrossSectionHandler->FindValue(Z, lowEdgeEnergy);
        ionloss   += e * cs * theAtomicNumDensityVector[iel];
      }	
      aVector->PutValue(i,ionloss);
    }
    theLossTable->insert(aVector);
  }
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














