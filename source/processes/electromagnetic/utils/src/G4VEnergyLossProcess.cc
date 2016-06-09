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
// $Id: G4VEnergyLossProcess.cc,v 1.2 2003/11/25 18:01:49 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEnergyLossProcess
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
// 15-02-03 Lambda table can be scaled (V.Ivanchenko)
// 17-02-03 Fix problem of store/restore tables (V.Ivanchenko)
// 18-02-03 Add control on CutCouple usage (V.Ivanchenko)
// 26-02-03 Simplify control on GenericIons (V.Ivanchenko)
// 06-03-03 Control on GenericIons using SubType + update verbose (V.Ivanchenko)
// 10-03-03 Add Ion registration (V.Ivanchenko)
// 22-03-03 Add Initialisation of cash (V.Ivanchenko)
// 26-03-03 Remove finalRange modification (V.Ivanchenko)
// 09-04-03 Fix problem of negative range limit for non integral (V.Ivanchenko)
// 26-04-03 Fix retrieve tables (V.Ivanchenko)
// 06-05-03 Set defalt finalRange = 1 mm (V.Ivanchenko)
// 12-05-03 Update range calculations + lowKinEnergy (V.Ivanchenko)
// 13-05-03 Add calculation of precise range (V.Ivanchenko)
// 23-05-03 Remove tracking cuts (V.Ivanchenko)
// 03-06-03 Fix initialisation problem for STD ionisation (V.Ivanchenko)
// 21-07-03 Add UpdateEmModel method (V.Ivanchenko)
// 03-11-03 Fix initialisation problem in RetrievePhysicsTable (V.Ivanchenko)
// 04-11-03 Add checks in RetrievePhysicsTable (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
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

#include "G4VEnergyLossProcess.hh"
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
#include "G4Proton.hh"
#include "G4VSubCutoffProcessor.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4GenericIon.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::G4VEnergyLossProcess(const G4String& name, G4ProcessType type):
                 G4VContinuousDiscreteProcess(name, type),
  nSCoffRegions(0),
  idxSCoffRegions(0),
  theDEDXTable(0),
  theRangeTable(0),
  theSecondaryRangeTable(0),
  theInverseRangeTable(0),
  theLambdaTable(0),
  theSubLambdaTable(0),
  theDEDXAtMaxEnergy(0),
  theRangeAtMaxEnergy(0),
  particle(0),
  baseParticle(0),
  secondaryParticle(0),
  currentCouple(0),
  nDEDXBins(90),
  nDEDXBinsForRange(70),
  nLambdaBins(90),
  faclow(1.5),
  linLossLimit(0.05),
  minSubRange(0.1),
  defaultRoverRange(0.2),
  defaultIntegralRange(1.0),
  lossFluctuationFlag(true),
  rndmStepFlag(false),
  hasRestProcess(true),
  tablesAreBuilt(false),
  integral(true),
  meanFreePath(true)
{

  minKinEnergy         = 0.1*keV;
  maxKinEnergy         = 100.0*GeV;
  maxKinEnergyForRange = 1.0*GeV;
  lowKinEnergy         = minKinEnergy*faclow;

  // default dRoverRange and finalRange
  SetStepFunction(defaultIntegralRange, 1.0*mm);
  SetVerboseLevel(0);

  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);
  scoffProcessors.clear();
  scoffRegions.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::~G4VEnergyLossProcess()
{
  Clear();

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
  delete modelManager;
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::Clear()
{
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::Clear() for " << GetProcessName() << G4endl;
  }
  if ( !baseParticle ) {
    if(theDEDXTable) theDEDXTable->clearAndDestroy();
    if(theRangeTable) theRangeTable->clearAndDestroy();
    if(theInverseRangeTable) theInverseRangeTable->clearAndDestroy();
    if(theLambdaTable) theLambdaTable->clearAndDestroy();
    if(theSubLambdaTable) theSubLambdaTable->clearAndDestroy();
    if(theDEDXAtMaxEnergy) delete [] theDEDXAtMaxEnergy;
    if(theRangeAtMaxEnergy) delete [] theRangeAtMaxEnergy;
  }

  theDEDXTable = 0;
  theRangeTable = 0;
  theInverseRangeTable = 0;
  theSecondaryRangeTable = 0;
  theLambdaTable = 0;
  theSubLambdaTable = 0;
  theDEDXAtMaxEnergy = 0;
  theRangeAtMaxEnergy = 0;
  tablesAreBuilt = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::Initialise()
{
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::Initialise() for "
           << GetProcessName() 
           << " for " << particle->GetParticleName()
           << G4endl;
  }

  Clear();

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

  const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if (nSCoffRegions) {
    const G4DataVector* theSubCuts = modelManager->SubCutoff();

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
    G4cout << "G4VEnergyLossProcess::Initialise() is done "
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

void G4VEnergyLossProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  currentCouple = 0;
  preStepLambda = 0.0;
  if(0 < verboseLevel) {
    G4cout << "========================================================" << G4endl;
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }

  if (part.GetParticleName() != "GenericIon" &&
      part.GetParticleType() == "nucleus" &&
      part.GetParticleSubType() == "generic")
  {
    (G4LossTableManager::Instance())->RegisterIon(&part, this);
    /*
    G4cout << part.GetProcessManager() << "  "
           << (G4GenericIon::GenericIon())->GetProcessManager()
           << G4endl;
    */
    return;
  }

  // Are particle defined?
  if( !particle ) {
    particle = &part;
    baseParticle = DefineBaseParticle(particle);
  }

  // Recalculation is needed because cuts were changed or recalculation is forced
  G4LossTableManager* lManager = G4LossTableManager::Instance();
  if ( lManager->IsRecalcNeeded(particle) ) {

    // It is responsability of the G4LossTables to build DEDX and range tables
    lManager->BuildPhysicsTable(particle);

    if(!baseParticle) PrintInfoDefinition();

    if(0 < verboseLevel) {
      G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() done for "
             << GetProcessName()
             << " and particle " << part.GetParticleName()
             << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetParticles(const G4ParticleDefinition* p1,
                                    const G4ParticleDefinition* p2)
{
  particle = p1;
  baseParticle = p2;
  G4bool yes = true;
  if(particle && baseParticle) {
    if((particle->GetPDGMass() < MeV && baseParticle->GetPDGMass() > MeV) ||
       (particle->GetPDGMass() > MeV && baseParticle->GetPDGMass() < MeV)) yes = false;

  }
  if(!yes) {
    G4cout << "Warning in G4VEnergyLossProcess::SetParticle: "
           << particle->GetParticleName()
           << " losses cannot be obtained from "
           << baseParticle->GetParticleName()
           << " losses" << G4endl;
  } else {
    DefineBaseParticle(p1);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::AddEmModel(G4int order, G4VEmModel* p, G4VEmFluctuationModel* fluc,
                                const G4Region* region)
{
  modelManager->AddEmModel(order, p, fluc, region);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::UpdateEmModel(const G4String& nam, G4double emin, G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::AddSubCutoffProcessor(G4VSubCutoffProcessor* p,
                                       const G4Region* r)
{
  if( !p ) {
    G4cout << "G4VEnergyLossProcess::AddSubCutoffProcessor WARNING: no SubCutoffProcessor defined." << G4endl;
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

G4PhysicsTable* G4VEnergyLossProcess::BuildDEDXTable()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable() for "
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
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable(): table is built for "
           << particle->GetParticleName()
           << G4endl;
    if(2 < verboseLevel) {
      G4cout << *theTable << G4endl;
    }
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildDEDXTableForPreciseRange()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTableForPreciseRange() for "
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
    G4PhysicsVector* aVector = DEDXPhysicsVectorForPreciseRange(couple);
    modelManager->FillDEDXVectorForPreciseRange(aVector, couple);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTableForPreciseRange(): table is built for " 
           << particle->GetParticleName()
           << G4endl;
    if(2 < verboseLevel) {
      G4cout << *theTable << G4endl;
    }
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaTable()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildLambdaTable() for process "
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
    G4cout << "Lambda table is built for "
           << particle->GetParticleName()
           << G4endl;
    if(2 < verboseLevel) {
      G4cout << *theTable << G4endl;
    }
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaSubTable()
{
  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildLambdaSubTable() for process "
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
    G4cout << "Table is built for "
           << particle->GetParticleName()
           << G4endl;
  }

  return theTable;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::AlongStepDoIt(const G4Track& track,
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
    const G4ParticleDefinition* d = track.GetDefinition();
    G4cout << "AlongStepDoIt for "
           << GetProcessName() << " and particle "
           << d->GetParticleName()
           << "  eScaled(MeV)= " << preStepScaledEnergy/MeV
           << "  slim(mm)= " << fRange/mm
           << "  s(mm)= " << length/mm
           << "  q^2= " << chargeSqRatio
           << " md= " << d->GetPDGMass()
           << G4endl;
  }
  */
    // low energy deposit case
  if (length >= fRange) {
    eloss = preStepKinEnergy;

  } else if(preStepScaledEnergy <= lowKinEnergy) {

    G4double x = 1.0 - length/fRange;
    eloss = preStepKinEnergy*(1.0 - x*x);

  // Short step
  } else if( length <= linLossLimit * fRange ) {
    eloss = (((*theDEDXTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))*length*chargeSqRatio;

  // Long step
  } else {
    G4double x = (fRange-length)/reduceFactor;
    G4double postStepScaledEnergy = ((*theInverseRangeTable)[currentMaterialIndex])->
             GetValue(x, b);
    eloss = (preStepScaledEnergy - postStepScaledEnergy)/massRatio;

    if (eloss <= 0.0) {
      eloss = (((*theDEDXTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))*length*chargeSqRatio;
    }

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

  }

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4double tmax = MaxSecondaryEnergy(dynParticle);
  tmax = std::min(tmax,(*theCuts)[currentMaterialIndex]);

  /*
  G4double eloss0 = eloss;
  if(-1 < verboseLevel) {
    //G4cout << *theDEDXTable << G4endl;
    G4cout << "eloss(MeV)= " << eloss/MeV
           << " eloss0(MeV)= " << (((*theDEDXTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))*length*chargeSqRatio
           << " r0(mm)= " << (((*theRangeTable)[currentMaterialIndex])->
               GetValue(preStepScaledEnergy, b))
           << " tmax= " << tmax
           << " e-eloss= " << preStepKinEnergy-eloss
      //   << " mat= " << currentMaterial
      //   << " matIdx= " << currentMaterialIndex
      //   << " preCouple= " << (step.GetPreStepPoint())->GetMaterialCutsCouple()
      //   << " postCouple= " << (step.GetPostStepPoint())->GetMaterialCutsCouple()
      //   << " rangeTable= " << theRangeTable
           << G4endl;
  }
  */

  // Sample fluctuations
  if (lossFluctuationFlag && eloss < preStepKinEnergy && eloss > 0.0) {

    eloss = modelManager->SampleFluctuations(currentMaterial, dynParticle,
                                       tmax, length, eloss, preStepScaledEnergy,
				       currentMaterialIndex);
  }

  /*
  if(-1 < verboseLevel) {
    G4cout << "eloss(MeV)= " << eloss/MeV
           << " fluc= " << (eloss-eloss0)/MeV
           << " currentChargeSquare= " << chargeSquare
           << " massRatio= " << massRatio
           << G4endl;
  }
  */

  // Subcutoff and/or deexcitation
  std::vector<G4Track*>* newp =
           SecondariesAlongStep(step, tmax, eloss, preStepScaledEnergy);

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
        if (pd != G4Gamma::Gamma() && pd != G4Electron::Electron() ) e += pd->GetPDGMass();

        preStepKinEnergy -= e;
        aParticleChange.AddSecondary(t);
      }
    }
    delete newp;
  }

  preStepKinEnergy -= eloss;

  /*
  if(-1 < verboseLevel) {
    G4cout << "eloss(MeV)= " << eloss/MeV
           << " preStepKinEnergy= " << preStepKinEnergy
           << " lossFlag= " << lossFluctuationFlag
           << G4endl;
  }
  */

  if (preStepKinEnergy <= 0.0) {

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

G4VParticleChange* G4VEnergyLossProcess::PostStepDoIt(const G4Track& track,
                                                  const G4Step& step)
{
  aParticleChange.Initialize(track);
  G4double finalT = track.GetKineticEnergy();
  G4double postStepScaledEnergy = finalT*massRatio;

  // Integral approach
  if (integral) {
    G4bool b;
    G4double postStepLambda = chargeSqRatio*
      (((*theLambdaTable)[currentMaterialIndex])->GetValue(postStepScaledEnergy,b));

    if(preStepLambda*G4UniformRand() > postStepLambda)
      return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
  }

  G4VEmModel* currentModel = SelectModel(postStepScaledEnergy);
  G4double tcut = (*theCuts)[currentMaterialIndex];
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();
  G4double tmax = currentModel->MaxSecondaryEnergy(dynParticle);

  /*
  if(0 < verboseLevel) {
    const G4ParticleDefinition* pd = dynParticle->GetDefinition();
    G4cout << "G4VEnergyLossProcess::PostStepDoIt: Sample secondary; E= " << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit(pd)
           << ", " <<  currentModel->HighEnergyLimit(pd) << ")"
           << G4endl;
  }
  */

  if (tcut < tmax)
    SecondariesPostStep(currentModel,currentCouple,dynParticle,tcut,finalT);

  if (finalT <= 0.0) {
    aParticleChange.SetEnergyChange(0.0);

    if (hasRestProcess) aParticleChange.SetStatusChange(fStopButAlive);
    else                aParticleChange.SetStatusChange(fStopAndKill);

    return pParticleChange;
  }

  aParticleChange.SetEnergyChange(finalT);

  return G4VContinuousDiscreteProcess::PostStepDoIt(track,step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::PrintInfoDefinition()
{
  G4cout << G4endl << GetProcessName() << ":   tables are built for  " 
         << particle->GetParticleName()
         << G4endl
         << "      dE/dx and range tables from "
	 << G4BestUnit(minKinEnergy,"Energy")
         << " to " << G4BestUnit(maxKinEnergy,"Energy")
         << " in " << nDEDXBins << " bins." << G4endl
         << "      Lambda tables from threshold to "
         << G4BestUnit(maxKinEnergy,"Energy")
         << " in " << nLambdaBins << " bins."
         << G4endl;
  if(theRangeTable) {
    G4cout << "      Step function: finalRange(mm)= " << finalRange/mm
           << ", dRoverRange= " << dRoverRange 
           << ", integral: " << integral
           << G4endl;
  }
  /*
      G4cout << "DEDXTable address= " << theDEDXTable << G4endl;
      if(theDEDXTable) G4cout << (*theDEDXTable) << G4endl;
      G4cout << "RangeTable address= " << theRangeTable << G4endl;
      if(theRangeTable) G4cout << (*theRangeTable) << G4endl;
      G4cout << "InverseRangeTable address= " << theInverseRangeTable << G4endl;
      if(theInverseRangeTable) G4cout << (*theInverseRangeTable) << G4endl;
  */
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
      G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
      G4cout << "SubLambdaTable address= " << theSubLambdaTable << G4endl;
      if(theSubLambdaTable) G4cout << (*theSubLambdaTable) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXTable(G4PhysicsTable* p)
{
  if(theDEDXTable) theDEDXTable->clearAndDestroy();
  theDEDXTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetRangeTable(G4PhysicsTable* p)
{
  if(theRangeTable) theRangeTable->clearAndDestroy();
  theRangeTable = p;
  size_t n = p->length();
  if(0 < n) {
    G4PhysicsVector* pv = (*p)[0];
    G4bool b;
    size_t nbins = pv->GetVectorLength();
    highKinEnergyForRange = pv->GetLowEdgeEnergy(nbins);
    theDEDXAtMaxEnergy = new G4double [n];
    theRangeAtMaxEnergy = new G4double [n];

    for (size_t i=0; i<n; i++) {
      pv = (*p)[i];
      G4double r2 = pv->GetValue(highKinEnergyForRange, b);
      G4double dedx = ((*theDEDXTable)[i])->GetValue(highKinEnergyForRange,b);
      theDEDXAtMaxEnergy[i] = dedx;
      theRangeAtMaxEnergy[i] = r2;
      //G4cout << "i= " << i << " e2= " << e2 << " r2= " << r2 
      //<< " dedx= " << dedx << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSecondaryRangeTable(G4PhysicsTable* p)
{
  theSecondaryRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetInverseRangeTable(G4PhysicsTable* p)
{
  if(theInverseRangeTable) theInverseRangeTable->clearAndDestroy();
  theInverseRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaTable(G4PhysicsTable* p)
{
  if(theLambdaTable) theLambdaTable->clearAndDestroy();
  theLambdaTable = p;
  tablesAreBuilt = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSubLambdaTable(G4PhysicsTable* p)
{
  if(theSubLambdaTable) theSubLambdaTable->clearAndDestroy();
  theSubLambdaTable = p;
  if (nSCoffRegions) {
    for (G4int i=0; i<nSCoffRegions; i++) {
      scoffProcessors[i]->SetLambdaSubTable(theSubLambdaTable);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::DEDXPhysicsVector(const G4MaterialCutsCouple* couple)
{
  G4int nbins = 3;
  if( couple->IsUsed() ) nbins = nDEDXBins;
  //G4double emax = maxKinEnergy*exp( log(maxKinEnergy/minKinEnergy) / ((G4double)(nbins-1)) );
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::DEDXPhysicsVectorForPreciseRange(
                             const G4MaterialCutsCouple* couple)
{
  G4int nbins = 3;
  if( couple->IsUsed() ) nbins = nDEDXBinsForRange;
  //G4double emax = maxKinEnergy*exp( log(maxKinEnergy/minKinEnergy) / ((G4double)(nbins-1)) );
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergyForRange, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::LambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  G4double cut  = (*theCuts)[couple->GetIndex()];
  G4int nbins = 3;
  if( couple->IsUsed() ) nbins = nLambdaBins;
  G4double tmin = std::max(MinPrimaryEnergy(particle, couple->GetMaterial(), cut),
                               minKinEnergy);
  if(tmin >= maxKinEnergy) tmin = 0.5*maxKinEnergy;
  //  G4double xmax = maxKinEnergy*exp(log(maxKinEnergy/tmin)/((G4double)(nbins-1)) );
  G4PhysicsVector* v = new G4PhysicsLogVector(tmin, maxKinEnergy, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEnergyLossProcess::SubLambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  return LambdaPhysicsVector(couple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MicroscopicCrossSection(G4double kineticEnergy,
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

G4double G4VEnergyLossProcess::MeanFreePath(const G4Track& track,
                                              G4double s,
                                              G4ForceCondition* cond)
{
  return GetMeanFreePath(track, s, cond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::ContinuousStepLimit(const G4Track& track,
                                               G4double x, G4double y, G4double& z)
{
  return GetContinuousStepLimit(track, x, y, z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetStepLimits(G4double v1, G4double v2)
{
  dRoverRange = v1;
  finalRange = v2;
  if (dRoverRange > 1.0) dRoverRange = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetIntegral(G4bool val)
{
  if(integral != val) {
    if(val) dRoverRange = defaultIntegralRange;
    else    dRoverRange = defaultRoverRange;
  }
  integral = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetStepFunction(G4double v1, G4double v2)
{
  dRoverRange = v1;
  finalRange = v2;
  if (dRoverRange > 0.999) dRoverRange = 1.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
  baseParticle = DefineBaseParticle(particle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetBaseParticle(const G4ParticleDefinition* p)
{
  baseParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::StorePhysicsTable(G4ParticleDefinition* part,
			 	     const G4String& directory,
				           G4bool ascii)
{
  G4bool res = true;
  if ( baseParticle ) return res;
  G4bool yes = true;

  if ( theDEDXTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
    yes = theDEDXTable->StorePhysicsTable(name,ascii);
    if( !yes ) res = false;
  }

  if ( theRangeTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Range",ascii);
    yes = theRangeTable->StorePhysicsTable(name,ascii);
    if( !yes ) res = false;
  }

  if ( theInverseRangeTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
    yes = theInverseRangeTable->StorePhysicsTable(name,ascii);
    if( !yes ) res = false;
  }

  if ( theLambdaTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    yes = theLambdaTable->StorePhysicsTable(name,ascii);
    if( !yes ) res = false;
  }

  if ( theSubLambdaTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"SubLambda",ascii);
    yes = theSubLambdaTable->StorePhysicsTable(name,ascii);
    if( !yes ) res = false;
  }
  if ( res ) {
    G4cout << "Physics tables are stored for " << particle->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  } else {
    G4cout << "Fail to store Physics Tables for " << particle->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VEnergyLossProcess::RetrievePhysicsTable(G4ParticleDefinition* part,
			  	        const G4String& directory,
			  	              G4bool ascii)
{
  G4bool res = true;
  currentCouple = 0;
  preStepLambda = 0.0;
  if(0 < verboseLevel) {
    G4cout << "========================================================" << G4endl;
    G4cout << "G4VEnergyLossProcess::RetrievePhysicsTable() for "
           << part->GetParticleName() << " and process " << GetProcessName() 
           << "; tables_are_built= " << tablesAreBuilt
           << G4endl;
  }

  const G4String particleName = part->GetParticleName();
  if( !particle ) {
    particle = part;
    baseParticle = DefineBaseParticle(particle);
  }

  if(particleName != "GenericIon"  &&
     part->GetParticleType() == "nucleus"  &&
     part->GetParticleSubType() == "generic")
  {
    (G4LossTableManager::Instance())->RegisterIon(part, this);
    return res;
  }

  if(tablesAreBuilt) return res;
  Initialise();

  // Recalculation is needed because cuts were changed or recalculation is forced
  G4LossTableManager* lManager = G4LossTableManager::Instance();
  if ( lManager->IsRecalcNeeded(particle)) {

    G4bool yes = true;
    G4bool fpi = true;
    if ( !baseParticle ) {
      G4PhysicsTable* table;
      G4String filename;
      const G4ProductionCutsTable* theCoupleTable=
            G4ProductionCutsTable::GetProductionCutsTable();
      size_t numOfCouples = theCoupleTable->GetTableSize();

      filename = GetPhysicsTableFileName(part,directory,"DEDX",ascii);
      
      table = new G4PhysicsTable(numOfCouples);
      yes = table->ExistPhysicsTable(filename);
      if(yes) yes = table->RetrievePhysicsTable(filename,ascii);
      if(yes) {
        SetDEDXTable(table);
        if (-1 < verboseLevel) {
          G4cout << "DEDX table for " << particleName << " is retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        fpi = false;
        table->clearAndDestroy();
        if (0 < verboseLevel) {
          G4cout << "DEDX table for " << particleName << " in file <"
                 << filename << "> is not retrieved"
                 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"Range",ascii);
      table = new G4PhysicsTable(numOfCouples);
      yes = table->ExistPhysicsTable(filename);
      if(yes) yes = table->RetrievePhysicsTable(filename,ascii);
      if(yes) {
        SetRangeTable(table);
        if (-1 < verboseLevel) {
          G4cout << "Range table for " << particleName << " is retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        table->clearAndDestroy();
        if(fpi) {
          res = false;
	  G4cout << "Range table for " << particleName << " in file <"
		 << filename << "> is not retrieved"
		 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"InverseRange",ascii);
      table = new G4PhysicsTable(numOfCouples);
      yes = table->ExistPhysicsTable(filename);
      if(yes)  yes = table->RetrievePhysicsTable(filename,ascii);
      if(yes) {
        SetInverseRangeTable(table);
        if (-1 < verboseLevel) {
          G4cout << "InverseRange table for " << particleName << " is retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        table->clearAndDestroy();
        if(fpi) {
          res = false;
          G4cout << "InverseRange table for " << particleName << " in file <"
                 << filename << "> is not retrieved"
                 << G4endl;

        }
      }

      filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
      table = new G4PhysicsTable(numOfCouples);
      yes = table->ExistPhysicsTable(filename);
      if(yes) yes = table->RetrievePhysicsTable(filename,ascii);
      if(yes) {
        SetLambdaTable(table);
        if (-1 < verboseLevel) {
          G4cout << "Lambda table for " << particleName << " is retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        table->clearAndDestroy();
        if(fpi) {
          res = false;
          G4cout << "Lambda table for " << particleName << " in file <"
                 << filename << "> is not retrieved"
                 << G4endl;
        }
      }

      filename = GetPhysicsTableFileName(part,directory,"SubLambda",ascii);
      table = new G4PhysicsTable(numOfCouples);
      yes = table->ExistPhysicsTable(filename);
      if(yes) yes = table->RetrievePhysicsTable(filename,ascii);
      if(yes) {
        SetSubLambdaTable(table);
        if (-1 < verboseLevel) {
          G4cout << "SubLambda table for " << particleName << " is retrieved from <"
                 << filename << ">"
                 << G4endl;
        }
      } else {
        table->clearAndDestroy();
        if(nSCoffRegions) {
          res=false;
          G4cout << "SubLambda table for " << particleName << " in file <"
                 << filename << "> is not retrieved"
                 << G4endl;
	}
      }
      if(res) PrintInfoDefinition();
      else {
        G4cout << "### BuildPhysicsTable will be requested for " <<  GetProcessName() 
               << " for " << particleName << G4endl; 
      }
    }
    tablesAreBuilt = true;
  }

  lManager->RetrievePhysicsTables(particle, this);

  return res;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
