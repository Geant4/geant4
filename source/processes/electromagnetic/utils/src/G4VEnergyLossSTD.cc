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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEnergyLossSTD
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 13-11-02 Minor fix - use normalised direction (V.Ivanchenko)
// 04-12-02 Minor change in PostStepDoIt (V.Ivanchenko)
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 26-12-02 Secondary production moved to derived classes (V.Ivanchenko)
// 04-01-03 Fix problem of very small steps for ions (V.Ivanchenko)
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 24-01-03 Temporarily close a control on usage of couples (V.Ivanchenko)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 05-02-03 Fix compilation warnings (V.Ivanchenko)
// 06-02-03 Add control on tmax in PostStepDoIt (V.Ivanchenko)
// 13-02-03 SubCutoffProcessors defined for regions (V.Ivanchenko)
//
// Class Description:
//
// It is the unified energy loss process it calculates the continuous
// energy loss for charged particles using a set of Energy Loss
// models valid for different energy regions. There are a possibility
// to create and access to dE/dx and range tables, or to calculate
// that information on fly.
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEnergyLossSTD.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4VSubCutoffProcessor.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossSTD::G4VEnergyLossSTD(const G4String& name, G4ProcessType type):
                 G4VContinuousDiscreteProcess(name, type),
  nSCoffRegions(0),
  idxSCoffRegions(0),
  theDEDXTable(0),
  theRangeTable(0),
  theSecondaryRangeTable(0),
  theInverseRangeTable(0),
  theLambdaTable(0),
  particle(0),
  baseParticle(0),
  secondaryParticle(0),
  theGamma(G4Gamma::Gamma()),
  theElectron(G4Electron::Electron()),
  minKinEnergy(1.0*eV),
  maxKinEnergy(100.0*GeV),
  linLossLimit(0.05),
  minSubRange(0.1),
  lossFluctuationFlag(true),
  rndmStepFlag(false),
  hasRestProcess(true),
  tablesAreBuilt(false),
  integral(true)
{
  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);
  scoffProcessors.clear();
  scoffRegions.clear();

  // default dRoverRange and finalRange
  SetStepLimits(0.2, 200.0*micrometer);

  verboseLevel = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossSTD::~G4VEnergyLossSTD()
{
  Clear();
  if(theLambdaTable) theLambdaTable->clearAndDestroy();
  delete modelManager;

  (G4LossTableManager::Instance())->Clear();

  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      if (scoffProcessors[i]) {
	for (G4int j=i+1; j<nSCoffRegions; j++) {
	  if(scoffProcessors[i] == scoffProcessors[j]) scoffProcessors[j] = 0;
	}
        delete scoffProcessors[i];
      }
    }
    scoffProcessors.clear();
    scoffRegions.clear();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::Clear()
{
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::Clear() for " << GetProcessName() << G4endl;
  }
  theDEDXTable = 0;
  theRangeTable = 0;
  theInverseRangeTable = 0;
  theSecondaryRangeTable = 0;
  modelManager->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::Initialise()
{

  Clear();

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::Initialise() for "
           << GetProcessName() << G4endl;
  }

  G4double initialCharge = particle->GetPDGCharge();
  G4double initialMass  = particle->GetPDGMass();
  chargeSquare = initialCharge*initialCharge/(eplus*eplus);
  chargeSqRatio = 1.0;
  massRatio = 1.0;
  reduceFactor = 1.0;

  if(particle->GetProcessManager()->GetAtRestProcessVector()->size())
               hasRestProcess = true;
  else         hasRestProcess = false;

  if (baseParticle) {
    massRatio = (baseParticle->GetPDGMass())/initialMass;
    G4double q = initialCharge/baseParticle->GetPDGCharge();
    chargeSqRatio = q*q;
    reduceFactor = 1.0/(chargeSqRatio*massRatio);
  }

  theCuts = modelManager->Initialise(particle, secondaryParticle, minSubRange, verboseLevel);
  
  // Sub Cutoff Regime

  idxSCoffRegions.clear();

  if (nSCoffRegions) {
    const G4DataVector* theSubCuts = modelManager->SubCutoff();
    const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();

    for (G4int i=0; i<nSCoffRegions; i++) {
      scoffProcessors[i]->Initialise(particle, secondaryParticle, theCuts, theSubCuts);
    }
    for (size_t j=0; j<numOfCouples; j++) {

      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(j);
      const G4ProductionCuts* pcuts = couple->GetProductionCuts();
      G4int reg = nSCoffRegions;
      do {reg--;} while (reg && pcuts != (scoffRegions[reg]->GetProductionCuts()));
      idxSCoffRegions.push_back(reg);
    }
  }
  if (0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::Initialise() is done "
           << " chargeSqRatio= " << chargeSqRatio
           << " massRatio= " << massRatio
           << " reduceFactor= " << reduceFactor << G4endl;
    if (nSCoffRegions) {
      G4cout << " SubCutoff Regime is ON for regions: " << G4endl;
      for (G4int i=0; i<nSCoffRegions; i++) {
        const G4Region* r = scoffRegions[i];
	G4cout << "           " << r->GetName() << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::BuildPhysicsTable(const G4ParticleDefinition& part)
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildPhysicsTable() for "
           << GetProcessName() << G4endl;
  }

  tablesAreBuilt = false;

  // Are particle defined?
  if(!particle || particle->GetParticleType() == "nucleus") {
    particle = &part;
  }
  if(particle != &part) {
    G4cout << "G4VEnergyLossSTD::BuildPhysicsTable: "
           << particle->GetParticleName()
           << " != " << part.GetParticleName()
           << G4endl;
    G4Exception("G4VEnergyLossSTD::BuildPhysicsTable with wrong G4ParticleDefinition");
  }

  if(!baseParticle) baseParticle = DefineBaseParticle(particle);

  if(particle->GetParticleType() == "nucleus" &&
     particle->GetParticleName() != "GenericIon" &&
     theLambdaTable) return;

  G4bool cutsWasModified = false;
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  for (size_t j=0; j<numOfCouples; j++){
    if (theCoupleTable->GetMaterialCutsCouple(j)->IsRecalcNeeded()) {
      cutsWasModified = true;
      break;
    }
  }
  if( !cutsWasModified ) return;

  // It is responsability of the G4LossTables to build DEDX and range tables
  (G4LossTableManager::Instance())->BuildDEDXTable(particle);


  //!!!! To be checked because finalRange is changed here !!!!
  for (size_t i=0; i<numOfCouples; i++) {
    G4double lengthCut = theCoupleTable->GetMaterialCutsCouple(i)
                       ->GetProductionCuts()->GetProductionCut(1);

    if (finalRange > lengthCut) finalRange = lengthCut;
  }

  if(theLambdaTable) theLambdaTable->clearAndDestroy();

  theLambdaTable = BuildLambdaTable();

  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      scoffProcessors[i]->SetLambdaSubTable(BuildLambdaSubTable());
    }
  }
  tablesAreBuilt = true;

  if(!baseParticle) PrintInfoDefinition();

  if(0 < verboseLevel) {
    G4cout << "Tables are built for " << particle->GetParticleName()
           << " IntegralFlag= " <<  integral
           << G4endl;
    if(2 < verboseLevel) {
      G4cout << "DEDXTable address= " << theDEDXTable << G4endl;
      if(theDEDXTable) G4cout << (*theDEDXTable) << G4endl;
      G4cout << "RangeTable address= " << theRangeTable << G4endl;
      if(theRangeTable) G4cout << (*theRangeTable) << G4endl;
      G4cout << "InverseRangeTable address= " << theInverseRangeTable << G4endl;
      if(theInverseRangeTable) G4cout << (*theInverseRangeTable) << G4endl;
    }
  }

  if(2 < verboseLevel) {
    G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
    if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossSTD::StorePhysicsTable(G4ParticleDefinition* part,
			 	     const G4String& directory,
				           G4bool ascii)
{
  G4String filename;

  if (!baseParticle) {

    G4String f1 = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
    G4String f2 = GetPhysicsTableFileName(part,directory,"Range",ascii);
    G4String f3 = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);

    if(!((G4LossTableManager::Instance())->StorePhysicsTable(this,f1,f2,f3,ascii)))
      {
        return false;
      }
  }

  if (theLambdaTable) {
    filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    if ( !theLambdaTable->StorePhysicsTable(filename, ascii) ){
      G4cout << "Fatal error theLambdaTable->StorePhysicsTable in <"
             << filename << ">"
             << G4endl;
      return false;
    }
  }

  G4PhysicsTable*  theLambdaSubTable = 0;
  if(nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      if(scoffProcessors[i]) {
        theLambdaSubTable = scoffProcessors[i]->LambdaSubTable();

        if (theLambdaSubTable) {
          filename = GetPhysicsTableFileName(part,directory,
	               "LambdaSub"+scoffProcessors[i]->GetName(),ascii);
          if ( !theLambdaSubTable->StorePhysicsTable(filename, ascii) ){
            G4cout << "Fatal error theLambdaSubTable->StorePhysicsTable in <"
                   << filename << ">"
                   << G4endl;
            return false;
	  }
        }
      }
    }
  }
  G4cout << GetProcessName() << " for " << particle->GetParticleName()
         << ": Success to store the PhysicsTables in "
         << directory << G4endl;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VEnergyLossSTD::RetrievePhysicsTable(G4ParticleDefinition* part,
				       const G4String& directory,
				             G4bool ascii)
{

  // Are particle defined?
  if(!particle) {
    particle = part;
  }
  if(particle != part) {
    G4cout << "G4VEnergyLossSTD::RetrievePhysicsTable: "
           << particle->GetParticleName()
           << " != " << part->GetParticleName()
           << G4endl;
    G4Exception("G4VEnergyLossSTD::RetrievePhysicsTable with wrong G4ParticleDefinition");
  }

  Initialise();

  G4PhysicsTable* theTable;
  G4String filename;
  G4String particleName = particle->GetParticleName();
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  filename = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
  theTable = new G4PhysicsTable(numOfCouples);
  if (theTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "DEDX table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;
  } else {
    G4cout << "DEDX table for " << particleName << " in file <"
           << filename << "> is not exist"
           << G4endl;
  }
  theDEDXTable = theTable;

  filename = GetPhysicsTableFileName(part,directory,"Range",ascii);
  theTable = new G4PhysicsTable(numOfCouples);
  if (theTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "Range table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;
  } else {
    G4cout << "Range table for " << particleName << " in file <"
           << filename << "> is not exist"
           << G4endl;
  }
  theRangeTable = theTable;

  filename = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
  theTable = new G4PhysicsTable(numOfCouples);
  if (theTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "Inverse table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;
    return false;
  }
  theInverseRangeTable = theTable;

  filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
  theLambdaTable = new G4PhysicsTable(numOfCouples);
  if (theLambdaTable->RetrievePhysicsTable(filename, ascii) ){
    G4cout << "Lambda table for " << particleName << " is retrieved from <"
           << filename << ">"
           << G4endl;
  } else {
    G4cout << "Lambda table for " << particleName << " in file <"
           << filename << "> is not exist"
           << G4endl;
  }

  G4PhysicsTable*  theLambdaSubTable = 0;
  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      if(scoffProcessors[i]) {

        filename = GetPhysicsTableFileName(part,directory,
                    "LambdaSub"+scoffProcessors[i]->GetName(),ascii);
        theLambdaSubTable = new G4PhysicsTable(numOfCouples);
        if (theLambdaSubTable->RetrievePhysicsTable(filename, ascii) ){
          G4cout << "LambdaSub table for " << particleName << " is retrieved from <"
                 << filename << ">"
                 << G4endl;

          scoffProcessors[i]->SetLambdaSubTable(theLambdaSubTable);

        } else {
          G4cout << "LambdaSub table for " << particleName << " in file <"
                 << filename << "> is not exist"
                 << G4endl;
        }
      }
    }
  }

  // G4LossTableManager manages DEDX and range tables

  (G4LossTableManager::Instance())->RetrieveDEDXTable(particle, this);

  G4cout << GetProcessName() << " for " << particleName
         << ": end of retrieving PhysicsTables from directory <"
         << directory << ">" << G4endl;

  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetParticles(const G4ParticleDefinition* p1,
                                    const G4ParticleDefinition* p2,
                                    const G4ParticleDefinition* p3)
{
  particle = p1;
  baseParticle = p2;
  secondaryParticle = p3;
  G4bool yes = false;
  if(particle) yes = particle->GetPDGCharge();
  if(baseParticle) {
    if(!baseParticle->GetPDGCharge()) yes = false;
    else {
      if((particle->GetPDGMass() < MeV &&
          baseParticle->GetPDGMass() > MeV) ||
         (particle->GetPDGMass() > MeV &&
          baseParticle->GetPDGMass() < MeV)) yes = false;
    }
  }
  if(!yes) {
    G4cout << "Warning in G4VEnergyLossSTD::SetParticle: "
           << particle->GetParticleName()
           << " losses cannot be obtained from "
           << baseParticle->GetParticleName()
           << " losses" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::AddEmModel(G4int order, G4VEmModel* p, G4VEmFluctuationModel* fluc,
                                const G4Region* region)
{
  modelManager->AddEmModel(order, p, fluc, region);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::AddSubCutoffProcessor(G4VSubCutoffProcessor* p,
                                       const G4Region* r)
{
  if( !p ) {
    G4cout << "G4VEnergyLossSTD::AddSubCutoffProcessor WARNING: no SubCutoffProcessor defined." << G4endl;
    return;
  }
  G4RegionStore* regionStore = G4RegionStore::GetInstance();
  if (!r) r = regionStore->GetRegion("DefaultRegionForTheWorld", false);
  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      if (r == scoffRegions[i]) {
        if ( scoffProcessors[i] ) delete scoffProcessors[i];
	scoffProcessors[i] = p;
        return;
      }
    }
  }
  scoffProcessors.push_back(p);
  scoffRegions.push_back(r);
  nSCoffRegions++;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossSTD::BuildDEDXTable()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildDEDXTable() for "
           << GetProcessName()
           << " and particle " << particle->GetParticleName()
           << G4endl;
  }

  // vectors to provide continues dE/dx
  G4DataVector factor;
  G4DataVector dedxLow;
  G4DataVector dedxHigh;

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  G4PhysicsTable* theTable = new G4PhysicsTable(numOfCouples);

  if(0 < verboseLevel) {
    G4cout << numOfCouples << " materials"
           << " minKinEnergy= " << minKinEnergy
           << " maxKinEnergy= " << maxKinEnergy
           << G4endl;
  }

  for(size_t i=0; i<numOfCouples; i++) {

    // create physics vector and fill it
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
    G4PhysicsVector* aVector = DEDXPhysicsVector(couple);
    modelManager->FillDEDXVector(aVector, couple);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildDEDXTable(): table is built" << G4endl;
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossSTD::BuildLambdaTable()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildLambdaTable() for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName()
           << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  G4PhysicsTable* theTable = new G4PhysicsTable(numOfCouples);

  for(size_t i=0; i<numOfCouples; i++) {

    // create physics vector and fill it
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
    G4PhysicsVector* aVector = LambdaPhysicsVector(couple);
    modelManager->FillLambdaVector(aVector, couple);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(0 < verboseLevel) {
    G4cout << "Lambda table is built" << G4endl;
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossSTD::BuildLambdaSubTable()
{
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildLambdaSubTable() for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName() << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  G4PhysicsTable* theTable = new G4PhysicsTable(numOfCouples);

  for(size_t i=0; i<numOfCouples; i++) {

    // create physics vector and fill it
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
    G4PhysicsVector* aVector = SubLambdaPhysicsVector(couple);
    modelManager->FillSubLambdaVector(aVector, couple);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(0 < verboseLevel) {
    G4cout << "Table is built" << G4endl;
  }

  return theTable;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossSTD::AlongStepDoIt(const G4Track& track,
                                                   const G4Step& step)
{
  aParticleChange.Initialize(track);
  // The process has range table - calculate energy loss
  if(!theRangeTable) return &aParticleChange;

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();
  G4double eloss  = 0.0;
  G4bool b;

  /*
  if(-1 < verboseLevel) {
    G4cout << "AlongStepDoIt for "
           << GetProcessName() << " and particle "
           << particle->GetParticleName()
           << "  eScaled(MeV)= " << preStepKinEnergy/MeV
           << "  slim(mm)= " << fRange/mm
           << "  s(mm)= " << length/mm
           << "  q^2= " << chargeSqRatio
           << G4endl;
  }
  */

  static const G4double faclow = 1.5;

    // low energy deposit case
  if (length >= fRange) {
    eloss = preStepKinEnergy;

  } else if(preStepScaledEnergy < faclow*minKinEnergy) {

    eloss = preStepKinEnergy*sqrt(length/fRange);

  // Short step
  } else if( length <= linLossLimit * fRange ) {
    eloss = (((*theDEDXTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))*length*chargeSqRatio;

  // Long step
  } else {
    G4double x = (fRange-length)/reduceFactor;
    G4PhysicsVector* v = (*theInverseRangeTable)[currentMaterialIndex];
    G4double postStepScaledEnergy = v->GetValue(x, b);
    eloss = (preStepScaledEnergy - postStepScaledEnergy)/massRatio;

    /*    
    if(-1 < verboseLevel) {
      G4cout << "fRange(mm)= " << fRange/mm
             << " xPost(mm)= " << x/mm
             << " ePre(MeV)= " << preStepScaledEnergy/MeV
             << " ePost(MeV)= " << postStepScaledEnergy/MeV
             << " eloss(MeV)= " << eloss/MeV
             << " eloss0(MeV)= " << (((*theDEDXTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))*length*chargeSqRatio
             << G4endl;
    }
    */

    if(eloss <= 0.0) eloss = preStepKinEnergy*sqrt(length/fRange);

  }

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4double tmax = MaxSecondaryEnergy(dynParticle);
  tmax = G4std::min(tmax,(*theCuts)[currentMaterialIndex]);

  // Sample fluctuations
  if (lossFluctuationFlag && eloss + minKinEnergy < preStepKinEnergy) {

    eloss = modelManager->SampleFluctuations(currentMaterial, dynParticle,
                                       tmax, length, eloss, preStepKinEnergy, 
				       currentMaterialIndex);
  }

  /*
  if(-1 < verboseLevel) {
    G4cout << "eloss(MeV)= " << eloss/MeV
           << " currentChargeSquare= " << chargeSquare
           << G4endl;
  }
  */  

  // Subcutoff and/or deexcitation
  G4std::vector<G4Track*>* newp =
           SecondariesAlongStep(step, tmax, eloss, preStepKinEnergy);

  if(newp) {

    G4int n = newp->size();
    if(n > 0) {
      aParticleChange.SetNumberOfSecondaries(n);
      G4Track* t;
      G4double e;
      for (G4int i=0; i<n; i++) {
        t = (*newp)[i];
        e = t->GetKineticEnergy();
        const G4ParticleDefinition* pd = t->GetDefinition();
        if (pd != theGamma && pd != theElectron ) e += pd->GetPDGMass();

        preStepKinEnergy -= e;
        aParticleChange.AddSecondary(t);
      }
    }
    delete newp;
  }

  /*
  if(-1 < verboseLevel) {
    G4cout << "eloss(MeV)= " << eloss/MeV
           << " preStepKinEnergy= " << preStepKinEnergy
           << " lossFlag= " << lossFluctuationFlag
           << G4endl;
  }
  */


  preStepKinEnergy -= eloss;

  if (preStepKinEnergy < minKinEnergy) {

    eloss += preStepKinEnergy;
    preStepKinEnergy = 0.0;

    if (hasRestProcess) aParticleChange.SetStatusChange(fStopButAlive);
    else                aParticleChange.SetStatusChange(fStopAndKill);
  }

  aParticleChange.SetEnergyChange(preStepKinEnergy);
  aParticleChange.SetLocalEnergyDeposit(eloss);

  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossSTD::PostStepDoIt(const G4Track& track,
                                                  const G4Step& step)
{
  aParticleChange.Initialize(track);
  G4double finalT = track.GetKineticEnergy();

  // Integral approach
  if(integral) {
    G4bool b;
    G4double postStepLambda =
      (((*theLambdaTable)[currentMaterialIndex])->GetValue(finalT,b))*chargeSqRatio;

    if(preStepLambda*G4UniformRand() > postStepLambda)
      return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
  }

  G4VEmModel* currentModel = SelectModel(finalT);
  G4double tcut = (*theCuts)[currentMaterialIndex];
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4double tmax = currentModel->MaxSecondaryEnergy(dynParticle);

  /*
  if(0 < verboseLevel) {
    const G4ParticleDefinition* pd = dynParticle->GetDefinition();
    G4cout << "G4VEnergyLossSTD::PostStepDoIt: Sample secondary; E= " << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit(pd)
           << ", " <<  currentModel->HighEnergyLimit(pd) << ")"
           << G4endl;
  }
  */

  if (tcut < tmax)
    SecondariesPostStep(currentModel,currentCouple,dynParticle,tcut,finalT);

  if (finalT < minKinEnergy) {
    aParticleChange.SetLocalEnergyDeposit(finalT);
    aParticleChange.SetEnergyChange(0.0);

    if (hasRestProcess) aParticleChange.SetStatusChange(fStopButAlive);
    else                aParticleChange.SetStatusChange(fStopAndKill);

    return pParticleChange;
  }

  aParticleChange.SetEnergyChange(finalT);

  return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::PrintInfoDefinition() const
{
  G4cout << G4endl << GetProcessName() << ":  " << G4endl
         << "      dE/dx and range tables from "
	 << G4BestUnit(MinKinEnergy(),"Energy")
         << " to " << G4BestUnit(MaxKinEnergy(),"Energy")
         << " in " << DEDXBinning() << " bins." << G4endl
         << "      Lambda tables from threshold to "
         << G4BestUnit(MaxKinEnergy(),"Energy")
         << " in " << LambdaBinning() << " bins."
         << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetDEDXTable(G4PhysicsTable* p)
{
  if(theDEDXTable) delete theDEDXTable;
  theDEDXTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetRangeTable(G4PhysicsTable* p)
{
  if(theRangeTable) delete theRangeTable;
  theRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetSecondaryRangeTable(G4PhysicsTable* p)
{
  if(theSecondaryRangeTable) delete theSecondaryRangeTable;
  theSecondaryRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossSTD::SetInverseRangeTable(G4PhysicsTable* p)
{
  if(theInverseRangeTable) delete theInverseRangeTable;
  theInverseRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossSTD::DEDXPhysicsVector(const G4MaterialCutsCouple* couple)
{
  //G4int nbins = 2;
  G4int nbins =  nDEDXBins;
  //if( couple->IsUsed() ) nbins = nDEDXBins;
  //  G4double xmax = maxKinEnergy*exp( log(maxKinEnergy/minKinEnergy) / ((G4double)(nbins-1)) );
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossSTD::LambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  G4double cut  = (*theCuts)[couple->GetIndex()];
  //G4int nbins = 2;
  G4int nbins = nLambdaBins;
  //if( couple->IsUsed() ) nbins = nLambdaBins;
  G4double tmin = G4std::max(MinPrimaryEnergy(particle, couple->GetMaterial(), cut),
                               minKinEnergy);
  if(tmin >= maxKinEnergy) tmin = 0.5*maxKinEnergy;
  //  G4double xmax = maxKinEnergy*exp(log(maxKinEnergy/tmin)/((G4double)(nbins-1)) );
  G4PhysicsVector* v = new G4PhysicsLogVector(tmin, maxKinEnergy, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossSTD::SubLambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  return LambdaPhysicsVector(couple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossSTD::MicroscopicCrossSection(G4double kineticEnergy,
                                             const G4MaterialCutsCouple* couple)
{
  // Cross section per atom is calculated
  DefineMaterial(couple);
  G4double cross = 0.0;
  G4bool b;
  if(theLambdaTable) {
    cross = (((*theLambdaTable)[currentMaterialIndex])->
                           GetValue(kineticEnergy, b));

    cross /= currentMaterial->GetTotNbOfAtomsPerVolume();
  }

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossSTD::MeanFreePath(const G4Track& track,
                                              G4double s,
                                              G4ForceCondition* cond)
{
  return GetMeanFreePath(track, s, cond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossSTD::ContinuousStepLimit(const G4Track& track,
                                               G4double x, G4double y, G4double& z)
{
  return GetContinuousStepLimit(track, x, y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

