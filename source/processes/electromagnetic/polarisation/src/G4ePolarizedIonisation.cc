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
// $Id: G4ePolarizedIonisation.cc,v 1.8 2010-06-16 11:20:54 schaelic Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4ePolarizedIonisation
//
// Author:        A.Schaelicke on base of Vladimir Ivanchenko code
//
// Creation date: 10.11.2005
//
// Modifications:
//
// 10-11-05, include polarization description (A.Schaelicke)
// , create asymmetry table and determine interactionlength 
// , update polarized differential cross section 
//
// 20-08-06, modified interface (A.Schaelicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
//
// Class Description:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ePolarizedIonisation.hh"
#include "G4Electron.hh"
#include "G4UniversalFluctuation.hh"
#include "G4BohrFluctuations.hh"
#include "G4UnitsTable.hh"

#include "G4PolarizedMollerBhabhaModel.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PolarizationManager.hh"
#include "G4PolarizationHelper.hh"
#include "G4StokesVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ePolarizedIonisation::G4ePolarizedIonisation(const G4String& name)
  : G4VEnergyLossProcess(name),
    theElectron(G4Electron::Electron()),
    isElectron(true),
    isInitialised(false),
    theAsymmetryTable(NULL),
    theTransverseAsymmetryTable(NULL)
{
  verboseLevel=0;
  //  SetDEDXBinning(120);
  //  SetLambdaBinning(120);
  //  numBinAsymmetryTable=78;

  //  SetMinKinEnergy(0.1*keV);
  //  SetMaxKinEnergy(100.0*TeV);
  //  PrintInfoDefinition();
  SetProcessSubType(fIonisation);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4ePolarizedIonisation::~G4ePolarizedIonisation()
{
  if (theAsymmetryTable) {
    theAsymmetryTable->clearAndDestroy();
    delete theAsymmetryTable;
  }
  if (theTransverseAsymmetryTable) {
    theTransverseAsymmetryTable->clearAndDestroy();
    delete theTransverseAsymmetryTable;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePolarizedIonisation::InitialiseEnergyLossProcess(
		    const G4ParticleDefinition* part,
		    const G4ParticleDefinition* /*part2*/)
{
  if(!isInitialised) {

    if(part == G4Positron::Positron()) isElectron = false;
    SetSecondaryParticle(theElectron);



    flucModel = new G4UniversalFluctuation();
    //flucModel = new G4BohrFluctuations();

    //    G4VEmModel* em = new G4MollerBhabhaModel();
    emModel = new G4PolarizedMollerBhabhaModel;
    emModel->SetLowEnergyLimit(MinKinEnergy());
    emModel->SetHighEnergyLimit(MaxKinEnergy());
    AddEmModel(1, emModel, flucModel);

    isInitialised = true;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4ePolarizedIonisation::PrintInfo()
{
  G4cout << "      Delta cross sections from Moller+Bhabha, "
         << "good description from 1 KeV to 100 GeV."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ePolarizedIonisation::GetMeanFreePath(const G4Track& track,
                                              G4double s,
                                              G4ForceCondition* cond)
{
  // *** get unploarised mean free path from lambda table ***
  G4double mfp = G4VEnergyLossProcess::GetMeanFreePath(track, s, cond);


  // *** get asymmetry, if target is polarized ***
  G4VPhysicalVolume*  aPVolume  = track.GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();

  G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();
  const G4bool volumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
  const G4StokesVector ePolarization = track.GetPolarization();

  if (mfp != DBL_MAX &&  volumeIsPolarized && !ePolarization.IsZero()) {
    const G4DynamicParticle* aDynamicElectron = track.GetDynamicParticle();
    const G4double eEnergy = aDynamicElectron->GetKineticEnergy();
    const G4ParticleMomentum eDirection0 = aDynamicElectron->GetMomentumDirection();

    G4StokesVector volumePolarization = polarizationManger->GetVolumePolarization(aLVolume);

    G4bool isOutRange;
    size_t idx = CurrentMaterialCutsCoupleIndex();
    G4double lAsymmetry = (*theAsymmetryTable)(idx)->
                                  GetValue(eEnergy, isOutRange);
    G4double tAsymmetry = (*theTransverseAsymmetryTable)(idx)->
                                  GetValue(eEnergy, isOutRange);

    // calculate longitudinal spin component
    G4double polZZ = ePolarization.z()*
			volumePolarization*eDirection0;
    // calculate transvers spin components
    G4double polXX = ePolarization.x()*
			volumePolarization*G4PolarizationHelper::GetParticleFrameX(eDirection0);
    G4double polYY = ePolarization.y()*
			volumePolarization*G4PolarizationHelper::GetParticleFrameY(eDirection0);


    G4double impact = 1. + polZZ*lAsymmetry + (polXX + polYY)*tAsymmetry;
    // determine polarization dependent mean free path
    mfp /= impact;
    if (mfp <=0.) {
     G4cout <<"PV impact ( "<<polXX<<" , "<<polYY<<" , "<<polZZ<<" )"<<G4endl;
     G4cout << " impact on MFP is "<< impact <<G4endl;
     G4cout<<" lAsymmetry= "<<lAsymmetry<<" ("<<std::fabs(lAsymmetry)-1.<<")\n";
     G4cout<<" tAsymmetry= "<<tAsymmetry<<" ("<<std::fabs(tAsymmetry)-1.<<")\n";
    }
  }

  return mfp;
}

G4double G4ePolarizedIonisation::PostStepGetPhysicalInteractionLength(const G4Track& track,
                                              G4double s,
                                              G4ForceCondition* cond)
{
  // *** get unploarised mean free path from lambda table ***
  G4double mfp = G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength(track, s, cond);


  // *** get asymmetry, if target is polarized ***
  G4VPhysicalVolume*  aPVolume  = track.GetVolume();
  G4LogicalVolume*    aLVolume  = aPVolume->GetLogicalVolume();

  G4PolarizationManager * polarizationManger = G4PolarizationManager::GetInstance();
  const G4bool volumeIsPolarized = polarizationManger->IsPolarized(aLVolume);
  const G4StokesVector ePolarization = track.GetPolarization();

  if (mfp != DBL_MAX &&  volumeIsPolarized && !ePolarization.IsZero()) {
    const G4DynamicParticle* aDynamicElectron = track.GetDynamicParticle();
    const G4double eEnergy = aDynamicElectron->GetKineticEnergy();
    const G4ParticleMomentum eDirection0 = aDynamicElectron->GetMomentumDirection();

    G4StokesVector volumePolarization = polarizationManger->GetVolumePolarization(aLVolume);

    G4bool isOutRange;
    size_t idx = CurrentMaterialCutsCoupleIndex();
    G4double lAsymmetry = (*theAsymmetryTable)(idx)->
                                  GetValue(eEnergy, isOutRange);
    G4double tAsymmetry = (*theTransverseAsymmetryTable)(idx)->
                                  GetValue(eEnergy, isOutRange);

    // calculate longitudinal spin component
    G4double polZZ = ePolarization.z()*
			volumePolarization*eDirection0;
    // calculate transvers spin components
    G4double polXX = ePolarization.x()*
			volumePolarization*G4PolarizationHelper::GetParticleFrameX(eDirection0);
    G4double polYY = ePolarization.y()*
			volumePolarization*G4PolarizationHelper::GetParticleFrameY(eDirection0);


    G4double impact = 1. + polZZ*lAsymmetry + (polXX + polYY)*tAsymmetry;
    // determine polarization dependent mean free path
    mfp /= impact;
    if (mfp <=0.) {
     G4cout <<"PV impact ( "<<polXX<<" , "<<polYY<<" , "<<polZZ<<" )"<<G4endl;
     G4cout << " impact on MFP is "<< impact <<G4endl;
     G4cout<<" lAsymmetry= "<<lAsymmetry<<" ("<<std::fabs(lAsymmetry)-1.<<")\n";
     G4cout<<" tAsymmetry= "<<tAsymmetry<<" ("<<std::fabs(tAsymmetry)-1.<<")\n";
    }
  }

  return mfp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void G4ePolarizedIonisation::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  // *** build DEDX and (unpolarized) cross section tables
  G4VEnergyLossProcess::BuildPhysicsTable(part);
  //  G4PhysicsTable* pt =
  //  BuildDEDXTable();


  // *** build asymmetry-table
  if (theAsymmetryTable) {
    theAsymmetryTable->clearAndDestroy(); delete theAsymmetryTable;}
  if (theTransverseAsymmetryTable) {
    theTransverseAsymmetryTable->clearAndDestroy(); delete theTransverseAsymmetryTable;}

  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  theAsymmetryTable = new G4PhysicsTable(numOfCouples);
  theTransverseAsymmetryTable = new G4PhysicsTable(numOfCouples);

  for (size_t j=0 ; j < numOfCouples; j++ ) {
    // get cut value
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);

    G4double cut = (*theCoupleTable->GetEnergyCutsVector(1))[j];

    //create physics vectors then fill it (same parameters as lambda vector)
    G4PhysicsVector * ptrVectorA = LambdaPhysicsVector(couple,cut);
    G4PhysicsVector * ptrVectorB = LambdaPhysicsVector(couple,cut);
    size_t nBins = ptrVectorA->GetVectorLength();

    for (size_t i = 0 ; i < nBins ; i++ ) {
      G4double lowEdgeEnergy = ptrVectorA->GetLowEdgeEnergy(i);
      G4double tasm=0.;
      G4double asym = ComputeAsymmetry(lowEdgeEnergy, couple, part, cut, tasm);
      ptrVectorA->PutValue(i,asym);
      ptrVectorB->PutValue(i,tasm);
    }
    theAsymmetryTable->insertAt( j , ptrVectorA ) ;
    theTransverseAsymmetryTable->insertAt( j , ptrVectorB ) ;
  }

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4ePolarizedIonisation::ComputeAsymmetry(G4double energy,
					 const G4MaterialCutsCouple* couple,
					       const G4ParticleDefinition& aParticle,
					       G4double cut,
					       G4double & tAsymmetry)
{
  G4double lAsymmetry = 0.0;
  	   tAsymmetry = 0.0;
  if (isElectron) {lAsymmetry = tAsymmetry = -1.0;}

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
  if (std::fabs(lAsymmetry)>1.) {
    G4cout<<" energy="<<energy<<"\n";
    G4cout<<"WARNING lAsymmetry= "<<lAsymmetry<<" ("<<std::fabs(lAsymmetry)-1.<<")\n";
  }
  if (std::fabs(tAsymmetry)>1.) {
    G4cout<<" energy="<<energy<<"\n";
    G4cout<<"WARNING tAsymmetry= "<<tAsymmetry<<" ("<<std::fabs(tAsymmetry)-1.<<")\n";
  }
//   else {
//     G4cout<<"        tAsymmetry= "<<tAsymmetry<<" ("<<std::fabs(tAsymmetry)-1.<<")\n";
//   }
  return lAsymmetry;
}


