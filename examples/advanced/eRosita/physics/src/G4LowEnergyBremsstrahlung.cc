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
// 
// --------------------------------------------------------------
//
// File name:     G4LowEnergyBremsstrahlung
//
// Author:        Alessandra Forti, Vladimir Ivanchenko
//
// Creation date: March 1999
//
// Modifications:
// 18.04.2000 V.L. 
//  - First implementation of continuous energy loss.
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
// 29.09.2001 V.Ivanchenko: major revision based on design iteration
// 10.10.2001 MGP Revision to improve code quality and consistency with design
// 18.10.2001 MGP Revision to improve code quality 
// 28.10.2001 VI  Update printout
// 29.11.2001 VI  New parametrisation
// 30.07.2002 VI  Fix in restricted energy loss
// 21.01.2003 VI  Cut per region
// 21.02.2003 V.Ivanchenko    Energy bins for spectrum are defined here
// 28.02.03 V.Ivanchenko    Filename is defined in the constructor
// 24.03.2003 P.Rodrigues Changes to accommodate new angular generators
// 20.05.2003 MGP  Removed memory leak related to angularDistribution
// 06.11.2003 MGP  Improved user interface to select angular distribution model
//
// --------------------------------------------------------------

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4RDeBremsstrahlungSpectrum.hh"
#include "G4RDBremsstrahlungCrossSectionHandler.hh"
#include "G4RDVBremAngularDistribution.hh"
#include "G4RDModifiedTsai.hh"
#include "G4RDGenerator2BS.hh"
#include "G4RDGenerator2BN.hh"
#include "G4RDVDataSetAlgorithm.hh"
#include "G4RDLogLogInterpolation.hh"
#include "G4RDVEMDataSet.hh"
#include "G4EnergyLossTables.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ProductionCutsTable.hh"


G4LowEnergyBremsstrahlung::G4LowEnergyBremsstrahlung(const G4String& nam)
  : G4eLowEnergyLoss(nam),
  crossSectionHandler(0),
  theMeanFreePath(0),
  energySpectrum(0)
{
  cutForPhotons = 0.;
  verboseLevel = 0;
  generatorName = "tsai";
  angularDistribution = new G4RDModifiedTsai("TsaiGenerator"); // default generator
//  angularDistribution->PrintGeneratorInformation();
  TsaiAngularDistribution = new G4RDModifiedTsai("TsaiGenerator");
}

/*
G4LowEnergyBremsstrahlung::G4LowEnergyBremsstrahlung(const G4String& nam, G4RDVBremAngularDistribution* distribution)
  : G4eLowEnergyLoss(nam),
    crossSectionHandler(0),
    theMeanFreePath(0),
    energySpectrum(0),
    angularDistribution(distribution)
{
  cutForPhotons = 0.;
  verboseLevel = 0;

  angularDistribution->PrintGeneratorInformation();

  TsaiAngularDistribution = new G4RDModifiedTsai("TsaiGenerator");
}
*/

G4LowEnergyBremsstrahlung::~G4LowEnergyBremsstrahlung()
{
  if(crossSectionHandler) delete crossSectionHandler;
  if(energySpectrum) delete energySpectrum;
  if(theMeanFreePath) delete theMeanFreePath;
  delete angularDistribution;
  delete TsaiAngularDistribution;
  energyBins.clear();
}


void G4LowEnergyBremsstrahlung::BuildPhysicsTable(const G4ParticleDefinition& aParticleType)
{
  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyBremsstrahlung::BuildPhysicsTable start"
           << G4endl;
      }

  cutForSecondaryPhotons.clear();

  // Create and fill BremsstrahlungParameters once
  if( energySpectrum != 0 ) delete energySpectrum;
  energyBins.clear();
  for(size_t i=0; i<15; i++) {
    G4double x = 0.1*((G4double)i);
    if(i == 0)  x = 0.01;
    if(i == 10) x = 0.95;
    if(i == 11) x = 0.97;
    if(i == 12) x = 0.99;
    if(i == 13) x = 0.995;
    if(i == 14) x = 1.0;
    energyBins.push_back(x);
  }
  const G4String dataName("/brem/br-sp.dat");
  energySpectrum = new G4RDeBremsstrahlungSpectrum(energyBins,dataName);

  if(verboseLevel > 0) {
    G4cout << "G4LowEnergyBremsstrahlungSpectrum is initialized"
           << G4endl;
      }

  // Create and fill G4RDCrossSectionHandler once

  if( crossSectionHandler != 0 ) delete crossSectionHandler;
  G4RDVDataSetAlgorithm* interpolation = new G4RDLogLogInterpolation();
  G4double lowKineticEnergy  = GetLowerBoundEloss();
  G4double highKineticEnergy = GetUpperBoundEloss();
  G4int    totBin = GetNbinEloss();
  crossSectionHandler = new G4RDBremsstrahlungCrossSectionHandler(energySpectrum, interpolation);
  crossSectionHandler->Initialise(0,lowKineticEnergy, highKineticEnergy, totBin);
  crossSectionHandler->LoadShellData("brem/br-cs-");

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
    G4cout << "G4LowEnergyBremsstrahlung::BuildPhysicsTable end"
           << G4endl;
      }
 
}


void G4LowEnergyBremsstrahlung::BuildLossTable(const G4ParticleDefinition& )
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

    // now comes the loop for the kinetic energy values
    for (size_t i = 0; i<totBin; i++) {

      G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(i);
      G4double ionloss = 0.;

      // loop for elements in the material
      for (size_t iel=0; iel<NumberOfElements; iel++ ) {
        G4int Z = (G4int)((*theElementVector)[iel]->GetZ());
        G4double e = energySpectrum->AverageEnergy(Z, 0.0, tCut, lowEdgeEnergy);
        G4double cs= crossSectionHandler->FindValue(Z, lowEdgeEnergy);
        ionloss   += e * cs  * theAtomicNumDensityVector[iel];
        if(verboseLevel > 1) {
          G4cout << "Z= " << Z
                 << "; tCut(keV)= " << tCut/keV
                 << "; E(keV)= " << lowEdgeEnergy/keV
                 << "; Eav(keV)= " << e/keV
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


G4VParticleChange* G4LowEnergyBremsstrahlung::PostStepDoIt(const G4Track& track,
							   const G4Step& step)
{
  aParticleChange.Initialize(track);

  const G4MaterialCutsCouple* couple = track.GetMaterialCutsCouple();
  G4double kineticEnergy = track.GetKineticEnergy();
  G4int index = couple->GetIndex();
  G4double tCut = cutForSecondaryPhotons[index];

  // Control limits
  if(tCut >= kineticEnergy)
     return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);

  G4int Z = crossSectionHandler->SelectRandomAtom(couple, kineticEnergy);

  G4double tGamma = energySpectrum->SampleEnergy(Z, tCut, kineticEnergy, kineticEnergy);
  G4double totalEnergy = kineticEnergy + electron_mass_c2;
  G4double finalEnergy = kineticEnergy - tGamma; // electron/positron final energy  
  G4double theta = 0;

  if((kineticEnergy < 1*MeV && kineticEnergy > 1*keV && generatorName == "2bn")){
      theta = angularDistribution->PolarAngle(kineticEnergy,finalEnergy,Z);
  }else{
      theta = TsaiAngularDistribution->PolarAngle(kineticEnergy,finalEnergy,Z);
  }

  G4double phi   = twopi * G4UniformRand();
  G4double dirZ  = std::cos(theta);
  G4double sinTheta  = std::sqrt(1. - dirZ*dirZ);
  G4double dirX  = sinTheta*std::cos(phi);
  G4double dirY  = sinTheta*std::sin(phi);

  G4ThreeVector gammaDirection (dirX, dirY, dirZ);
  G4ThreeVector electronDirection = track.GetMomentumDirection();

  //
  // Update the incident particle 
  //
  gammaDirection.rotateUz(electronDirection);   
    
  // Kinematic problem
  if (finalEnergy < 0.) {
    tGamma += finalEnergy;
    finalEnergy = 0.0;
  }

  G4double momentum = std::sqrt((totalEnergy + electron_mass_c2)*kineticEnergy);

  G4double finalX = momentum*electronDirection.x() - tGamma*gammaDirection.x();
  G4double finalY = momentum*electronDirection.y() - tGamma*gammaDirection.y();
  G4double finalZ = momentum*electronDirection.z() - tGamma*gammaDirection.z();
      
  aParticleChange.SetNumberOfSecondaries(1);
  G4double norm = 1./std::sqrt(finalX*finalX + finalY*finalY + finalZ*finalZ);
  aParticleChange.ProposeMomentumDirection(finalX*norm, finalY*norm, finalZ*norm);
  aParticleChange.ProposeEnergy( finalEnergy );

  // create G4DynamicParticle object for the gamma 
  G4DynamicParticle* aGamma= new G4DynamicParticle (G4Gamma::Gamma(),
						    gammaDirection, tGamma);
  aParticleChange.AddSecondary(aGamma);

  return G4VContinuousDiscreteProcess::PostStepDoIt(track, step);
}


void G4LowEnergyBremsstrahlung::PrintInfoDefinition()
{
  G4String comments = "Total cross sections from EEDL database.";
  comments += "\n      Gamma energy sampled from a parameterised formula.";
  comments += "\n      Implementation of the continuous dE/dx part.";  
  comments += "\n      At present it can be used for electrons ";
  comments += "in the energy range [250eV,100GeV].";
  comments += "\n      The process must work with G4LowEnergyIonisation.";
  
  G4cout << G4endl << GetProcessName() << ":  " << comments << G4endl;
}         

G4bool G4LowEnergyBremsstrahlung::IsApplicable(const G4ParticleDefinition& particle)
{
  return (  (&particle == G4Electron::Electron())  );
}


G4double G4LowEnergyBremsstrahlung::GetMeanFreePath(const G4Track& track,
						    G4double,
						    G4ForceCondition* cond)
{
  *cond = NotForced;
  G4int index = (track.GetMaterialCutsCouple())->GetIndex();
  const G4RDVEMDataSet* data = theMeanFreePath->GetComponent(index);
  G4double meanFreePath = data->FindValue(track.GetKineticEnergy());
  return meanFreePath;
}

void G4LowEnergyBremsstrahlung::SetCutForLowEnSecPhotons(G4double cut)
{
  cutForPhotons = cut;
}

void G4LowEnergyBremsstrahlung::SetAngularGenerator(G4RDVBremAngularDistribution* distribution)
{
  angularDistribution = distribution;
  angularDistribution->PrintGeneratorInformation();
}

void G4LowEnergyBremsstrahlung::SetAngularGenerator(const G4String& name)
{
  if (name == "tsai") 
    {
      delete angularDistribution;
      angularDistribution = new G4RDModifiedTsai("TsaiGenerator");
      generatorName = name;
    }
  else if (name == "2bn")
    {
      delete angularDistribution;
      angularDistribution = new G4RDGenerator2BN("2BNGenerator");
      generatorName = name;
    }
  else if (name == "2bs")
    {
       delete angularDistribution;
       angularDistribution = new G4RDGenerator2BS("2BSGenerator");
       generatorName = name;
    }
  else
    {
      G4Exception("G4LowEnergyBremsstrahlung::SetAngularGenerator()",
                  "InvalidSetup", FatalException, "Generator does not exist!");
    }

  angularDistribution->PrintGeneratorInformation();
}

