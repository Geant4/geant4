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
// $Id: G4eplusPolarizedAnnihilation.cc 105740 2017-08-16 13:05:44Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4eplusPolarizedAnnihilation
//
// Author:        A. Schaelicke on base of Vladimir Ivanchenko / Michel Maire code
//
// Creation date: 02.07.2006
//
// Modifications:
// 26-07-06 modified cross section  (P. Starovoitov)
// 21-08-06 interface updated   (A. Schaelicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
// 02-10-07, enable AtRest (V.Ivanchenko)
//
//
// Class Description:
//
// Polarized process of e+ annihilation into 2 gammas
//

//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4eplusPolarizedAnnihilation.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Gamma.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"

#include "G4PolarizedAnnihilationModel.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusPolarizedAnnihilation::G4eplusPolarizedAnnihilation(const G4String& name)
  : G4eplusAnnihilation(name), isInitialised(false),
    theAsymmetryTable(nullptr),
    theTransverseAsymmetryTable(nullptr)
{
  emModel = new G4PolarizedAnnihilationModel();
  SetEmModel(emModel); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4eplusPolarizedAnnihilation::~G4eplusPolarizedAnnihilation()
{
  CleanTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::CleanTables()
{
  if(theAsymmetryTable) {
    theAsymmetryTable->clearAndDestroy();
    delete theAsymmetryTable;
    theAsymmetryTable = nullptr;
  }
  if(theTransverseAsymmetryTable) {
    theTransverseAsymmetryTable->clearAndDestroy();
    delete theTransverseAsymmetryTable;
    theTransverseAsymmetryTable = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eplusPolarizedAnnihilation::GetMeanFreePath(const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition)
{
  G4double mfp = G4VEmProcess::GetMeanFreePath(track, previousStepSize, condition);

  if(theAsymmetryTable && theTransverseAsymmetryTable && mfp < DBL_MAX) {
    mfp *= ComputeSaturationFactor(track);
  }
  if (verboseLevel>=2) {
    G4cout << "G4eplusPolarizedAnnihilation::MeanFreePath:  " 
	   << mfp / mm << " mm " << G4endl;
  }
  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eplusPolarizedAnnihilation::PostStepGetPhysicalInteractionLength(
                              const G4Track& track,
                              G4double previousStepSize,
                              G4ForceCondition* condition)
{
  // save previous values
  G4double nLength = theNumberOfInteractionLengthLeft;
  G4double iLength = currentInteractionLength;

  // *** compute unpolarized step limit ***
  // this changes theNumberOfInteractionLengthLeft and currentInteractionLength
  G4double x = G4VEmProcess::PostStepGetPhysicalInteractionLength(track, 
								  previousStepSize, 
								  condition);
  G4double x0 = x;
  G4double satFact = 1.0;
 
  // *** add corrections on polarisation ***
  if(theAsymmetryTable && theTransverseAsymmetryTable && x < DBL_MAX) {
    satFact = ComputeSaturationFactor(track);
    G4double curLength = currentInteractionLength*satFact;
    G4double prvLength = iLength*satFact;
    if(nLength > 0.0) {
      theNumberOfInteractionLengthLeft = 
	std::max(nLength - previousStepSize/prvLength, 0.0);
    }
    x = theNumberOfInteractionLengthLeft * curLength;
  }
  if (verboseLevel>=2) {
    G4cout << "G4eplusPolarizedAnnihilation::PostStepGPIL: "
           << std::setprecision(8) << x/mm  << " mm;" << G4endl
           << "                         unpolarized value: "
           << std::setprecision(8) << x0/mm << " mm." << G4endl;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4eplusPolarizedAnnihilation::ComputeSaturationFactor(const G4Track& track)
{
  G4Material*         aMaterial = track.GetMaterial();
  G4VPhysicalVolume*  aPVolume  = track.GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();
    
  G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();
    
  const G4bool volumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
  G4StokesVector electronPolarization = polarizationManger->GetVolumePolarization(aLVolume);

  G4double factor = 1.0;

  if (volumeIsPolarized) {
     
    // *** get asymmetry, if target is polarized ***
    const G4DynamicParticle* aDynamicPositron = track.GetDynamicParticle();
    const G4double positronEnergy = aDynamicPositron->GetKineticEnergy();
    const G4StokesVector positronPolarization = track.GetPolarization();
    const G4ParticleMomentum positronDirection0 = aDynamicPositron->GetMomentumDirection();

    if (verboseLevel>=2) {
      G4cout << "G4eplusPolarizedAnnihilation::ComputeSaturationFactor: " << G4endl;      
      G4cout << " Mom " << positronDirection0  << G4endl;
      G4cout << " Polarization " << positronPolarization  << G4endl;
      G4cout << " MaterialPol. " << electronPolarization  << G4endl;
      G4cout << " Phys. Volume " << aPVolume->GetName() << G4endl;
      G4cout << " Log. Volume  " << aLVolume->GetName() << G4endl;
      G4cout << " Material     " << aMaterial          << G4endl;
    }
    
    size_t midx = CurrentMaterialCutsCoupleIndex();
    const G4PhysicsVector* aVector = nullptr;
    const G4PhysicsVector* bVector = nullptr;
    if(midx < theAsymmetryTable->size()) { 
      aVector = (*theAsymmetryTable)(midx);
    }
    if(midx < theTransverseAsymmetryTable->size()) { 
      bVector = (*theTransverseAsymmetryTable)(midx);
    }
    if(aVector && bVector) {
      G4double lAsymmetry = aVector->Value(positronEnergy);
      G4double tAsymmetry = bVector->Value(positronEnergy);
      G4double polZZ = positronPolarization.z()*
	(electronPolarization*positronDirection0);
      G4double polXX = positronPolarization.x()*
	(electronPolarization*G4PolarizationHelper::GetParticleFrameX(positronDirection0));
      G4double polYY = positronPolarization.y()*
	(electronPolarization*G4PolarizationHelper::GetParticleFrameY(positronDirection0));

      factor /= (1. + polZZ*lAsymmetry + (polXX + polYY)*tAsymmetry);

      if (verboseLevel>=2) {
	G4cout << " Asymmetry:     " << lAsymmetry << ", " << tAsymmetry  << G4endl;
	G4cout << " PolProduct:    " << polXX << ", " << polYY << ", " << polZZ << G4endl;
	G4cout << " Factor:        " << factor         << G4endl;
      }
    } else {
      G4ExceptionDescription ed;
      ed << "Problem with asymmetry tables: material index " << midx 
	 << " is out of range or tables are not filled";
      G4Exception("G4eplusPolarizedAnnihilation::ComputeSaturationFactor","em0048",
		  JustWarning, ed, "");
    }
  }
  return factor;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::BuildPhysicsTable(
	    const G4ParticleDefinition& part) 
{
  G4VEmProcess::BuildPhysicsTable(part);
  G4bool isMaster = true;
  const  G4eplusPolarizedAnnihilation* masterProcess = 
    static_cast<const G4eplusPolarizedAnnihilation*>(GetMasterProcess());
  if(masterProcess && masterProcess != this) { isMaster = false; }
  if(isMaster) { BuildAsymmetryTables(part); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::BuildAsymmetryTables(
            const G4ParticleDefinition& part)
{
  // cleanup old, initialise new table
  CleanTables();
  theAsymmetryTable = 
    G4PhysicsTableHelper::PreparePhysicsTable(theAsymmetryTable);
  theTransverseAsymmetryTable = 
    G4PhysicsTableHelper::PreparePhysicsTable(theTransverseAsymmetryTable);

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  //G4cout<<" annih-numOfCouples="<<numOfCouples<<"\n";
  for(size_t i=0; i<numOfCouples; ++i) {
    //G4cout<<"annih- "<<i<<"/"<<numOfCouples<<"\n";
    if (!theAsymmetryTable) break;
    //G4cout<<"annih- "<<theAsymmetryTable->GetFlag(i)<<"\n";
    if (theAsymmetryTable->GetFlag(i)) {
      //G4cout<<" building pol-annih ... \n";

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);

      // use same parameters as for lambda
      G4PhysicsVector* aVector = LambdaPhysicsVector(couple);
      G4PhysicsVector* tVector = LambdaPhysicsVector(couple);

      for (G4int j = 0 ; j < LambdaBinning() ; ++j ) {
	G4double lowEdgeEnergy = aVector->GetLowEdgeEnergy(j);
	G4double tasm=0.;
	G4double asym = ComputeAsymmetry(lowEdgeEnergy, couple, part, 0., tasm);
	aVector->PutValue(j,asym);
	tVector->PutValue(j,tasm);
      }
      G4PhysicsTableHelper::SetPhysicsVector(theAsymmetryTable, i, aVector);
      G4PhysicsTableHelper::SetPhysicsVector(theTransverseAsymmetryTable, i, tVector);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4eplusPolarizedAnnihilation::ComputeAsymmetry(G4double energy,
			    const G4MaterialCutsCouple* couple,
			    const G4ParticleDefinition& aParticle,
			    G4double cut,
			    G4double &tAsymmetry)
{
  G4double lAsymmetry = 0.0;
  tAsymmetry = 0.0;

  // calculate polarized cross section
  theTargetPolarization=G4ThreeVector(0.,0.,1.);
  emModel->SetTargetPolarization(theTargetPolarization);
  emModel->SetBeamPolarization(theTargetPolarization);
  G4double sigma2=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

  // calculate transversely polarized cross section
  theTargetPolarization=G4ThreeVector(1.,0.,0.);
  emModel->SetTargetPolarization(theTargetPolarization);
  emModel->SetBeamPolarization(theTargetPolarization);
  G4double sigma3=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

  // calculate unpolarized cross section
  theTargetPolarization=G4ThreeVector();
  emModel->SetTargetPolarization(theTargetPolarization);
  emModel->SetBeamPolarization(theTargetPolarization);
  G4double sigma0=emModel->CrossSection(couple,&aParticle,energy,cut,energy);

  // determine assymmetries
  if (sigma0>0.) {
    lAsymmetry=sigma2/sigma0-1.;
    tAsymmetry=sigma3/sigma0-1.;
  }
  return lAsymmetry;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4eplusPolarizedAnnihilation::PrintInfo()
{
  G4cout << "      Polarized model for annihilation into 2 photons"
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
