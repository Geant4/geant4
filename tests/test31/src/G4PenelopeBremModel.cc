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
// $Id: G4PenelopeBremModel.cc,v 1.1 2004-07-27 08:33:04 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------
//
// File name:     G4PenelopeBremModel
//
// Author:        Luciano Pandola
//
// Creation date: February 2003
//
// Modifications:
// 24.04.2003 V.Ivanchenko - Cut per region mfpt
// 20.05.2003 MGP          - Removed compilation warnings
//                           Restored NotForced in GetMeanFreePath
// 23.05.2003 L.Pandola    - Bug fixed. G4double cast to elements of
//                           ebins vector (is it necessary?)
// 23.05.2003 MGP          - Removed memory leak (fix in destructor)
//
//----------------------------------------------------------------

#include "G4PenelopeBremModel.hh"
#include "G4PenelopeBremsstrahlungContinuous.hh"
#include "G4eBremsstrahlungSpectrum.hh"
#include "G4BremsstrahlungCrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DataVector.hh"
#include "G4ProductionCutsTable.hh"
#include "Randomize.hh"


//
//
// Il load dei dati continui si puo' mettere nel BuildLossTable
//
//
// Anche la sezione d'urto totale dovrebbe essere convertita
// le tabelle di G4 mi paiono piu' complete pero' (le lascio!)
//

G4PenelopeBremModel::G4PenelopeBremModel(const G4ParticleDefinition* p, const G4String& nam)
  : G4VEmModel(nam),
  crossSectionHandler(0),
  theMeanFreePath(0),
  energySpectrum(0),
  particle(0),
  highKinEnergy(100.*TeV),
  lowKinEnergy(1.0*keV),
  minThreshold(0.1*keV),
  totBins(70)
{
  particle = p;
  verboseLevel = 0;
}


G4PenelopeBremModel::~G4PenelopeBremModel()
{
  delete crossSectionHandler;
  delete energySpectrum;
  delete theMeanFreePath;

  for (size_t m=0; m<materialAngularData.size(); m++)
    {
      G4AngularData materialData = *(materialAngularData[m]);
      for (size_t i=0; i<materialData.size(); i++)
      {
	delete materialData[i];
      }
      delete &materialData;
    }
}


void G4PenelopeBremModel::Initialise(const G4ParticleDefinition* aParticleType,
                                                           const G4DataVector& cutForSecondaryPhotons)
{
  if(verboseLevel > 0) {
    G4cout << "G4PenelopeBremModel::Initialise start"
           << G4endl;
  }

  LoadAngularData();
  // Create and fill BremModelParameters once
  if ( energySpectrum != 0 ) delete energySpectrum;
  //grid of reduced energy bins for photons

  G4DataVector eBins;

  //G4double value;
  eBins.push_back(1.0e-12);
  eBins.push_back(0.05);
  eBins.push_back(0.075);
  eBins.push_back(0.1);
  eBins.push_back(0.125);
  eBins.push_back(0.15);
  eBins.push_back(0.2);
  eBins.push_back(0.25);
  eBins.push_back(0.3);
  eBins.push_back(0.35);
  eBins.push_back(0.40);
  eBins.push_back(0.45);
  eBins.push_back(0.50);
  eBins.push_back(0.55);
  eBins.push_back(0.60);
  eBins.push_back(0.65);
  eBins.push_back(0.70);
  eBins.push_back(0.75);
  eBins.push_back(0.80);
  eBins.push_back(0.85);
  eBins.push_back(0.90);
  eBins.push_back(0.925);
  eBins.push_back(0.95);
  eBins.push_back(0.97);
  eBins.push_back(0.99);
  eBins.push_back(0.995);
  eBins.push_back(0.999);
  eBins.push_back(0.9995);
  eBins.push_back(0.9999);
  eBins.push_back(0.99995);
  eBins.push_back(0.99999);
  eBins.push_back(1.0);


  const G4String dataName("/penelope/br-sp-pen.dat");
  energySpectrum = new G4eBremsstrahlungSpectrum(eBins,dataName);

  //the shape of the energy spectrum for positron is the same used for the electrons,
  //as the differential cross section is scaled of a factor f(E,Z) which is independent
  //on the energy of the gamma

  if(verboseLevel > 0) {
    G4cout << "G4PenelopeBremModelSpectrum is initialized"
           << G4endl;
        }

  particle = aParticleType;
  // Create and fill G4CrossSectionHandler once

  if( crossSectionHandler != 0 ) delete crossSectionHandler;
  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();

  crossSectionHandler = new G4BremsstrahlungCrossSectionHandler(energySpectrum, interpolation);
  crossSectionHandler->Initialise(0,lowKinEnergy, highKinEnergy, totBins);
  if (aParticleType==G4Electron::Electron())
    {
      //G4cout << "Leggo le sezioni d'urto degli elettroni" << G4endl;
      crossSectionHandler->LoadShellData("brem/br-cs-");
    }
  else
    {
      //G4cout << "Leggo le sezioni d'urto dei positroni" << G4endl;
      crossSectionHandler->LoadShellData("penelope/br-cs-pos-"); //cross section for positrons
    }

  if (verboseLevel > 0) {
    G4cout << "### Cross section data: "
           << G4endl;
    crossSectionHandler->PrintData();
    G4cout << "Parameters: "
           << G4endl;
    energySpectrum->PrintData();
    }

  if( theMeanFreePath != 0 ) delete theMeanFreePath;
  theMeanFreePath = crossSectionHandler->
                    BuildMeanFreePathForMaterials(&cutForSecondaryPhotons);

  if(verboseLevel > 0) {
    G4cout << "G4PenelopeBremModel::Initialise end"
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremModel::HighEnergyLimit(const G4ParticleDefinition*)
{
  return highKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremModel::LowEnergyLimit(const G4ParticleDefinition*)
{
  return lowKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4PenelopeBremModel::MinEnergyCut(const G4ParticleDefinition*,
                                              const G4MaterialCutsCouple*)
{
  return minThreshold;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4PenelopeBremModel::IsInCharge(const G4ParticleDefinition* p)
{
  return (p == G4Electron::Electron() || p == G4Positron::Positron());
}


G4double G4PenelopeBremModel::ComputeDEDX(const G4Material* material,
                                             const G4ParticleDefinition* ,
                                                   G4double kineticEnergy,
                                                   G4double cutEnergy)
{
  // Build table for energy loss due to soft brems
  // the tables are built for *MATERIALS* binning is taken from LowEnergyLoss

    G4double tCut = std::min(highKinEnergy, cutEnergy);

    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector =
      material->GetAtomicNumDensityVector();
    if(verboseLevel > 1) {
      G4cout << "Energy loss for material " << material->GetName()
             << " tCut(keV)= " << cutEnergy/keV
             << G4endl;
      }

    G4double loss = 0.0;

    const G4String partName = particle->GetParticleName();
    // loop for elements in the material
    for (size_t iel=0; iel<NumberOfElements; iel++ ) {
      G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
      //costruttore del continuo per l'elemento Z
      G4PenelopeBremsstrahlungContinuous* ContLoss =
           new G4PenelopeBremsstrahlungContinuous(Z,tCut,kineticEnergy*0.5,kineticEnergy,partName);
      loss += ContLoss->CalculateStopping(kineticEnergy)  * theAtomicNumDensityVector[iel];
      delete ContLoss;
    }
    return loss;
}

G4DynamicParticle* G4PenelopeBremModel::SampleSecondary(
                                const G4MaterialCutsCouple* couple,
                                const G4DynamicParticle* dp,
                                      G4double tmin,
                                      G4double maxEnergy)
{
  const G4Material* material = couple->GetMaterial();
  G4double kineticEnergy = dp->GetKineticEnergy();
  G4int index = couple->GetIndex();
  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);
  G4double tGamma = energySpectrum->SampleEnergy(Z, tmin, maxEnergy, kineticEnergy);
  //bisogna recuperare il puntatore
  G4AngularData* elementData =  materialAngularData[index]; //punta al vettore che contiene i puntatori del materiale
  //look for the index of the selected atom
  const G4ElementVector* theElementVector = material->GetElementVector();
  size_t NumberOfElements = material->GetNumberOfElements();
  // loop for elements in the material
  size_t iel=0;
  G4int Z_try=0,indexEl=0;
  do {
    if (iel >= NumberOfElements ) {
      G4String excep = "Not found the angular data for material " + material->GetName();
      G4Exception(excep);
    }
    Z_try = (G4int)((*theElementVector)[(size_t) iel]->GetZ());
    indexEl = iel;
    iel++;
  }while(Z_try != Z);

  G4PenelopeBremsstrahlungAngular* finalAngularData = (*elementData)[indexEl];

  // Sample gamma angle (Z - axis along the parent particle).
  G4double dirZ = finalAngularData->ExtractCosTheta(kineticEnergy,tGamma);

  //std::ofstream fff("prova.dat",std::ios::app);
  //fff << dirZ << G4endl;
  //fff.close();

  G4double phi   = twopi * G4UniformRand();
  G4double sinTheta  = sqrt(1. - dirZ*dirZ);
  G4double dirX  = sinTheta*cos(phi);
  G4double dirY  = sinTheta*sin(phi);

  G4ThreeVector gammaDirection (dirX, dirY, dirZ);


  // create G4DynamicParticle object for the gamma
  G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
						    gammaDirection, tGamma);

  return aGamma;
}


G4double G4PenelopeBremModel::CrossSection(const G4Material* material,
                        const G4ParticleDefinition*,
                              G4double kineticEnergy,
                              G4double,
                              G4double)
{
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();

  G4int index = 0;
  for (G4int m=0; m<numOfCouples; m++) {
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
    if (material == couple->GetMaterial()) {
      index = m;
      break;
    }
  }

  const G4VEMDataSet* data = theMeanFreePath->GetComponent(index);
  G4double meanFreePath = data->FindValue(kineticEnergy);
  G4double cross = 0.0;
  if(meanFreePath<DBL_MAX) cross = 1.0/meanFreePath;
  return cross;
}

void G4PenelopeBremModel::LoadAngularData()
{
  materialAngularData.clear();

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  G4int numOfCouples = theCoupleTable->GetTableSize();

  for (G4int m=0; m<numOfCouples; m++) {
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(m);
    // get material parameters needed for the energy loss calculation
    const G4Material* material= couple->GetMaterial();
    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements();
    // loop for elements in the material
    G4AngularData* elementAngularData = new G4AngularData;
    for (size_t iel=0; iel<NumberOfElements; iel++ ) {
     G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
     G4PenelopeBremsstrahlungAngular* AngPointer = new G4PenelopeBremsstrahlungAngular(Z);
     elementAngularData->push_back(AngPointer);
    }
    materialAngularData.push_back(elementAngularData);
  }
}
