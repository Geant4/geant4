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
// $Id: G4PenelopeBremsstrahlung.cc,v 1.10 2003/06/16 17:00:16 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
// 
// --------------------------------------------------------------
//
// File name:     G4PenelopeBremsstrahlung
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

#include "G4PenelopeBremsstrahlung.hh"
#include "G4PenelopeBremsstrahlungContinuous.hh"
#include "G4eBremsstrahlungSpectrum.hh"
#include "G4BremsstrahlungCrossSectionHandler.hh"
#include "G4VDataSetAlgorithm.hh"
#include "G4LogLogInterpolation.hh"
#include "G4VEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4DataVector.hh"
#include "G4ProductionCutsTable.hh"


//
//
// Il load dei dati continui si puo' mettere nel BuildLossTable
//
//
// Anche la sezione d'urto totale dovrebbe essere convertita
// le tabelle di G4 mi paiono piu' complete pero' (le lascio!)
//

G4PenelopeBremsstrahlung::G4PenelopeBremsstrahlung(const G4String& nam)
  : G4eLowEnergyLoss(nam), 
  crossSectionHandler(0),
  theMeanFreePath(0),
  energySpectrum(0)
{
  // materialAngularData.clear();
  cutForPhotons = 0.;
  verboseLevel = 0;
  LoadAngularData();
}


G4PenelopeBremsstrahlung::~G4PenelopeBremsstrahlung()
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


void G4PenelopeBremsstrahlung::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4PenelopeBremsstrahlung::BuildPhysicsTable start"
           << G4endl;
  }

  cutForSecondaryPhotons.clear();

  // Create and fill BremsstrahlungParameters once
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
    G4cout << "G4PenelopeBremsstrahlungSpectrum is initialized"
           << G4endl;
        }

  // Create and fill G4CrossSectionHandler once

  if( crossSectionHandler != 0 ) delete crossSectionHandler;
  G4VDataSetAlgorithm* interpolation = new G4LogLogInterpolation();
  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  G4int    totBin = GetNbinEloss();
  crossSectionHandler = new G4BremsstrahlungCrossSectionHandler(energySpectrum, interpolation);
  crossSectionHandler->Initialise(0,lowKineticEnergy, highKineticEnergy, totBin);
  if (&aParticleType==G4Electron::Electron()) 
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
    G4cout << GetProcessName() 
           << " is created; Cross section data: " 
           << G4endl;
    crossSectionHandler->PrintData();
    G4cout << "Parameters: " 
           << G4endl;
    energySpectrum->PrintData();
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
  // Build mean free path data using cut values

  if( theMeanFreePath != 0 ) delete theMeanFreePath;
  theMeanFreePath = crossSectionHandler->
                    BuildMeanFreePathForMaterials(&cutForSecondaryPhotons);

  if(verboseLevel > 0) {
    G4cout << "The MeanFreePath table is built"
           << G4endl;
      }

  // Build common DEDX table for all ionisation processes
 
  BuildDEDXTable(aParticleType);

  if(verboseLevel > 0) {
    G4cout << "G4PenelopeBremsstrahlung::BuildPhysicsTable end"
           << G4endl;
  }
}


void G4PenelopeBremsstrahlung::BuildLossTable(const G4ParticleDefinition& aParticleType)
{
  // Build table for energy loss due to soft brems
  // the tables are built for *MATERIALS* binning is taken from LowEnergyLoss
  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  size_t totBin = GetNbinEloss();
 
  //  create table
  
  if (theLossTable) { 
      theLossTable->clearAndDestroy();
      delete theLossTable;
  }

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  theLossTable = new G4PhysicsTable(numOfCouples);
 
  
  // Clean up the vector of cuts
  cutForSecondaryPhotons.clear();

  // Loop for materials
  for (size_t j=0; j<numOfCouples; j++) {
    
    // create physics vector and fill it
    G4PhysicsLogVector* aVector = new G4PhysicsLogVector(lowKineticEnergy,
		        				 highKineticEnergy,
							 totBin);

    // get material parameters needed for the energy loss calculation
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
    const G4Material* material= couple->GetMaterial();
    // the cut cannot be below lowest limit
    G4double tCut = (*(theCoupleTable->GetEnergyCutsVector(0)))[j];
    tCut = std::min(highKineticEnergy, tCut);
    cutForSecondaryPhotons.push_back(tCut);
   
    const G4ElementVector* theElementVector = material->GetElementVector();
    size_t NumberOfElements = material->GetNumberOfElements() ;
    const G4double* theAtomicNumDensityVector = 
      material->GetAtomicNumDensityVector();
    if(verboseLevel > 1) {
      G4cout << "Energy loss for material # " << j
             << " tCut(keV)= " << tCut/keV
             << G4endl;
      }

    G4DataVector* ionloss = new G4DataVector();
    for (size_t i = 0; i<totBin; i++) {
      ionloss->push_back(0.0); //il vettore ha come elementi totBin zeri
    }
            
    const G4String partName = aParticleType.GetParticleName();
    // loop for elements in the material
    for (size_t iel=0; iel<NumberOfElements; iel++ ) {
      G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
      //costruttore del continuo per l'elemento Z
      G4PenelopeBremsstrahlungContinuous* ContLoss = new G4PenelopeBremsstrahlungContinuous(Z,tCut,GetLowerBoundEloss(),
											    GetUpperBoundEloss(),
											    partName);
      //l'inizializzazione va ripetuta per ogni Z!
      // now comes the loop for the kinetic energy values
      for (size_t k = 0; k<totBin; k++) {
	G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(k);
	//qui viene chiamato il metodo giusto della perdita continua
    	(*ionloss)[k] += ContLoss->CalculateStopping(lowEdgeEnergy)  * theAtomicNumDensityVector[iel];
	//if (Z == 29)
	//  G4cout << lowEdgeEnergy << " " << ContLoss->CalculateStopping(lowEdgeEnergy)  * theAtomicNumDensityVector[iel] << G4endl;
	//va chiamato una volta per ogni energia
	//per ogni energia k, somma su tutti gli elementi del materiale
      }
      delete ContLoss;
    }

    for (size_t ibin = 0; ibin<totBin; ibin++) {
      aVector->PutValue(ibin,(*ionloss)[ibin]); //un valore per ogni energia (somma sugli elementi del mate)
    }
    delete ionloss;
    theLossTable->insert(aVector);
  }
}


G4VParticleChange* G4PenelopeBremsstrahlung::PostStepDoIt(const G4Track& track,
							  const G4Step& step)
{
  aParticleChange.Initialize(track);

  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  const G4Material* material = couple->GetMaterial();
  G4double kineticEnergy = track.GetKineticEnergy();
  G4int index = couple->GetIndex();
  G4double tCut = cutForSecondaryPhotons[index];

  // Control limits
  if(tCut >= kineticEnergy)
     return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);

  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);
  G4double tGamma = energySpectrum->SampleEnergy(Z, tCut, kineticEnergy, kineticEnergy);
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

  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double phi   = twopi * G4UniformRand();
  G4double sinTheta  = sqrt(1. - dirZ*dirZ);
  G4double dirX  = sinTheta*cos(phi);
  G4double dirY  = sinTheta*sin(phi); 
    
  G4ThreeVector gammaDirection (dirX, dirY, dirZ);
  G4ThreeVector electronDirection = track.GetMomentumDirection();
    
  gammaDirection.rotateUz(electronDirection);   
  
  //
  // Update the incident particle 
  //
    
  G4double finalEnergy = kineticEnergy - tGamma;  
    
  // Kinematic problem
  if (finalEnergy < 0.) {
    tGamma += finalEnergy;
    finalEnergy = 0.0;
  }

  G4double momentum = sqrt((totalEnergy + electron_mass_c2)*kineticEnergy);

  G4double finalX = momentum*electronDirection.x() - tGamma*gammaDirection.x();
  G4double finalY = momentum*electronDirection.y() - tGamma*gammaDirection.y();
  G4double finalZ = momentum*electronDirection.z() - tGamma*gammaDirection.z();

  aParticleChange.SetNumberOfSecondaries(1);
  G4double norm = 1./sqrt(finalX*finalX + finalY*finalY + finalZ*finalZ); 
  aParticleChange.SetMomentumChange(finalX*norm, finalY*norm, finalZ*norm);
  aParticleChange.SetEnergyChange( finalEnergy );

  // create G4DynamicParticle object for the gamma 
  G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
						    gammaDirection, tGamma);
  aParticleChange.AddSecondary(aGamma); 

  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


void G4PenelopeBremsstrahlung::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database for electrons.";
  comments += "\n Total cross section for positrons calculated from the electrons";
  comments += "\n through an empirical scaling function.";
  comments += "\n      Gamma energy sampled from a data-driven histogram.";
  comments += "\n      Implementation of the continuous dE/dx part.";  
  comments += "\n      It can be used for electrons and positrons";
  comments += "in the energy range [250eV,100GeV].";
  comments += "\n      The process must work with G4LowEnergyIonisation.";

  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}

G4bool G4PenelopeBremsstrahlung::IsApplicable(const G4ParticleDefinition& particle)
{
  return (  (&particle == G4Electron::Electron()) || (&particle == G4Positron::Positron()) );
}


G4double G4PenelopeBremsstrahlung::GetMeanFreePath(const G4Track& track,
						   G4double, // previousStepSize
						    G4ForceCondition* cond)
{
  *cond = NotForced;
  G4int index = (track.GetMaterialCutsCouple())->GetIndex();
  const G4VEMDataSet* data = theMeanFreePath->GetComponent(index);
  G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
  return meanFreePath; 
} 

void G4PenelopeBremsstrahlung::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForPhotons = cut;
}

void G4PenelopeBremsstrahlung::LoadAngularData()
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  const size_t numOfMaterials = G4Material::GetNumberOfMaterials();
  materialAngularData.clear();
  
  for (size_t j=0; j<numOfMaterials; j++) {
    // get material parameters needed for the energy loss calculation
    const G4Material* material= (*theMaterialTable)[j];
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
