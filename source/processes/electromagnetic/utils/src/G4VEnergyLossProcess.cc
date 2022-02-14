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
// Modifications: Vladimir Ivanchenko
//
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4LossTableManager.hh"
#include "G4LossTableBuilder.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4EmParameters.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4PhysicsModelCatalog.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4SafetyHelper.hh"
#include "G4TransportationManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4VSubCutProducer.hh"
#include "G4EmBiasingManager.hh"
#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::G4VEnergyLossProcess(const G4String& name, 
                                           G4ProcessType type): 
  G4VContinuousDiscreteProcess(name, type)
{
  theParameters = G4EmParameters::Instance();
  SetVerboseLevel(1);

  // low energy limit
  lowestKinEnergy = theParameters->LowestElectronEnergy();

  // Size of tables
  minKinEnergy     = 0.1*CLHEP::keV;
  maxKinEnergy     = 100.0*CLHEP::TeV;
  maxKinEnergyCSDA = 1.0*CLHEP::GeV;
  nBins            = 84;
  nBinsCSDA        = 35;

  // default linear loss limit
  finalRange = 1.*CLHEP::mm;

  // particle types
  theElectron   = G4Electron::Electron();
  thePositron   = G4Positron::Positron();
  theGamma      = G4Gamma::Gamma();

  // run time objects
  pParticleChange = &fParticleChange;
  fParticleChange.SetSecondaryWeightByProcess(true);
  modelManager = new G4EmModelManager();
  safetyHelper = G4TransportationManager::GetTransportationManager()
    ->GetSafetyHelper();
  aGPILSelection = CandidateForSelection;

  // initialise model
  lManager = G4LossTableManager::Instance();
  lManager->Register(this);

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  theDensityFactor = bld->GetDensityFactors();
  theDensityIdx = bld->GetCoupleIndexes();

  scTracks.reserve(10);
  secParticles.reserve(12);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::~G4VEnergyLossProcess()
{
  /*
  G4cout << "** G4VEnergyLossProcess::~G4VEnergyLossProcess() for " 
         << GetProcessName() << " isMaster: " << isMaster
         << "  basePart: " << baseParticle 
         << G4endl;
  G4cout << " isIonisation " << isIonisation << "  " 
         << theDEDXTable << "  " <<  theIonisationTable << G4endl;
  */

  if (isMaster && nullptr == baseParticle) {
    if(nullptr != theDEDXTable) {

      //G4cout << " theIonisationTable " << theIonisationTable << G4endl;
      if(theIonisationTable == theDEDXTable) { theIonisationTable = nullptr; }
      //G4cout << " delete theDEDXTable " << theDEDXTable << G4endl;
      theDEDXTable->clearAndDestroy();
      delete theDEDXTable;
      theDEDXTable = nullptr;
    }
    //G4cout << " theIonisationTable " << theIonisationTable << G4endl;
    if(nullptr != theIonisationTable) {
      //G4cout << " delete theIonisationTable " << theIonisationTable << G4endl;
      theIonisationTable->clearAndDestroy();
      delete theIonisationTable;
      theIonisationTable = nullptr;
    }
    if(nullptr != theDEDXunRestrictedTable && isIonisation) {
      theDEDXunRestrictedTable->clearAndDestroy();
      delete theDEDXunRestrictedTable;
      theDEDXunRestrictedTable = nullptr;
    }
    if(nullptr != theCSDARangeTable && isIonisation) {
      theCSDARangeTable->clearAndDestroy();
      delete theCSDARangeTable;
      theCSDARangeTable = nullptr;
    }
    //G4cout << "delete RangeTable: " << theRangeTableForLoss << G4endl;
    if(nullptr != theRangeTableForLoss && isIonisation) {
      theRangeTableForLoss->clearAndDestroy();
      delete theRangeTableForLoss;
      theRangeTableForLoss = nullptr;
    }
    //G4cout << "delete InvRangeTable: " << theInverseRangeTable << G4endl;
    if(nullptr != theInverseRangeTable && isIonisation /*&& !isIon*/) {
      theInverseRangeTable->clearAndDestroy();
      delete theInverseRangeTable;
      theInverseRangeTable = nullptr;
    }
    //G4cout << "delete LambdaTable: " << theLambdaTable << G4endl;
    if(nullptr != theLambdaTable) {
      theLambdaTable->clearAndDestroy();
      delete theLambdaTable;
      theLambdaTable = nullptr;
    }
    if(nullptr != fXSpeaks) {
      for(auto const & v : *fXSpeaks) { delete v; }
      delete fXSpeaks;
      fXSpeaks = nullptr;
    }
  }
  secParticles.clear();
  delete modelManager;
  delete biasManager;
  delete scoffRegions;
  delete emModels;
  lManager->DeRegister(this);
  //G4cout << "** all removed" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MinPrimaryEnergy(const G4ParticleDefinition*, 
                                                const G4Material*, 
                                                G4double cut)
{
  return cut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::AddEmModel(G4int order, G4VEmModel* ptr,
                                      G4VEmFluctuationModel* fluc,
                                      const G4Region* region)
{
  if(nullptr == ptr) { return; }
  modelManager->AddEmModel(order, ptr, fluc, region);
  ptr->SetParticleChange(pParticleChange, fluc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetEmModel(G4VEmModel* ptr, G4int)
{
  if(nullptr == ptr) { return; }
  if(nullptr == emModels) { emModels = new std::vector<G4VEmModel*>; }
  if(!emModels->empty()) {
    for(auto & em : *emModels) { if(em == ptr) { return; } }
  }
  emModels->push_back(ptr);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDynamicMassCharge(G4double massratio,
                                                G4double charge2ratio)
{
  massRatio = massratio;
  logMassRatio = G4Log(massRatio);
  fFactor = charge2ratio*biasFactor;
  if(baseMat) { fFactor *= (*theDensityFactor)[currentCoupleIndex]; }
  chargeSqRatio = charge2ratio;
  reduceFactor  = 1.0/(fFactor*massRatio);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::PreparePhysicsTable for "
           << GetProcessName() << " for " << part.GetParticleName() 
           << "  " << this << G4endl;
  }
  isMaster = lManager->IsMaster();

  // Are particle defined?
  if(nullptr == particle) { particle = &part; }

  if(part.GetParticleType() == "nucleus") {

    G4String pname = part.GetParticleName();
    if(pname != "deuteron" && pname != "triton" &&
       pname != "alpha+"   && pname != "alpha") {

      if(nullptr == theGenericIon) {
        theGenericIon = 
          G4ParticleTable::GetParticleTable()->FindParticle("GenericIon");
      }
      isIon = true; 
      if(particle != theGenericIon) {
        G4ProcessManager* pm = theGenericIon->GetProcessManager();
        G4ProcessVector* v = pm->GetAlongStepProcessVector();
        size_t n = v->size();
        for(size_t j=0; j<n; ++j) {
          if((*v)[j] == this) {
            particle = theGenericIon;
            break;
          } 
        }
      }
    }
  }

  if( particle != &part ) {
    if(!isIon) {
      lManager->RegisterExtraParticle(&part, this);
    }
    if(1 < verboseLevel) {
      G4cout << "### G4VEnergyLossProcess::PreparePhysicsTable()"
             << " interrupted for "
             << part.GetParticleName() << "  isIon=" << isIon
             << " baseMat=" << baseMat 
             << "  particle " << particle << "  GenericIon " << theGenericIon 
             << G4endl;
    }
    return;
  }

  tablesAreBuilt = false;

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  lManager->PreparePhysicsTable(&part, this, isMaster);

  // Base particle and set of models can be defined here
  InitialiseEnergyLossProcess(particle, baseParticle);

  // parameters of the process
  if(!actLossFluc) { lossFluctuationFlag = theParameters->LossFluctuation(); }
  rndmStepFlag = theParameters->UseCutAsFinalRange();
  if(!actMinKinEnergy) { minKinEnergy = theParameters->MinKinEnergy(); }
  if(!actMaxKinEnergy) { maxKinEnergy = theParameters->MaxKinEnergy(); }
  if(!actBinning) { nBins = theParameters->NumberOfBins(); }
  maxKinEnergyCSDA = theParameters->MaxEnergyForCSDARange();
  nBinsCSDA = theParameters->NumberOfBinsPerDecade()
    *G4lrint(std::log10(maxKinEnergyCSDA/minKinEnergy));
  if(!actLinLossLimit) { linLossLimit = theParameters->LinearLossLimit(); }
  lambdaFactor    = theParameters->LambdaFactor();
  logLambdafactor = G4Log(lambdaFactor);
  if(isMaster) { SetVerboseLevel(theParameters->Verbose()); }
  else { SetVerboseLevel(theParameters->WorkerVerbose()); }

  theParameters->DefineRegParamForLoss(this);

  fRangeEnergy = fLambdaEnergy = 0.0;

  G4double initialCharge = particle->GetPDGCharge();
  G4double initialMass   = particle->GetPDGMass();

  theParameters->FillStepFunction(particle, this);

  // integral option may be disabled
  if(!theParameters->Integral()) { fXSType = fEmNoIntegral; }

  // parameters for scaling from the base particle
  if (nullptr != baseParticle) {
    massRatio    = (baseParticle->GetPDGMass())/initialMass;
    logMassRatio = G4Log(massRatio);
    G4double q = initialCharge/baseParticle->GetPDGCharge();
    chargeSqRatio = q*q;
    if(chargeSqRatio > 0.0) { reduceFactor = 1.0/(chargeSqRatio*massRatio); }
  }
  lowestKinEnergy = (initialMass < CLHEP::MeV) 
    ? theParameters->LowestElectronEnergy()
    : theParameters->LowestMuHadEnergy();

  // Tables preparation
  if (isMaster && nullptr == baseParticle) {

    if(nullptr != theDEDXTable && isIonisation) {
      if(nullptr != theIonisationTable && theDEDXTable != theIonisationTable) {
        theDEDXTable->clearAndDestroy();
        delete theDEDXTable;
        theDEDXTable = theIonisationTable;
      }   
    }
    
    theDEDXTable = G4PhysicsTableHelper::PreparePhysicsTable(theDEDXTable);
    bld->InitialiseBaseMaterials(theDEDXTable);

    if (theParameters->BuildCSDARange()) {
      theDEDXunRestrictedTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theDEDXunRestrictedTable);
      theCSDARangeTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theCSDARangeTable);
    }

    theLambdaTable = G4PhysicsTableHelper::PreparePhysicsTable(theLambdaTable);

    if(isIonisation) {
      theRangeTableForLoss = 
        G4PhysicsTableHelper::PreparePhysicsTable(theRangeTableForLoss);
      theInverseRangeTable = 
        G4PhysicsTableHelper::PreparePhysicsTable(theInverseRangeTable);  
    }

    if(fXSType == fEmTwoPeaks) {
      const G4ProductionCutsTable* theCoupleTable=
	G4ProductionCutsTable::GetProductionCutsTable();
      size_t n = theCoupleTable->GetTableSize();
      if(nullptr == fXSpeaks) { 
	fXSpeaks = new std::vector<G4TwoPeaksXS*>;
      }
      fXSpeaks->resize(n, nullptr);
    }
  }
  /*
  G4cout << "** G4VEnergyLossProcess::PreparePhysicsTable() for " 
         << GetProcessName() << " and " << particle->GetParticleName()
         << " isMaster: " << isMaster << " isIonisation: " 
         << isIonisation << G4endl;
  G4cout << " theDEDX: " << theDEDXTable 
         << " theRange: " << theRangeTableForLoss
         << " theInverse: " << theInverseRangeTable
         << " theLambda: " << theLambdaTable << G4endl;
  */
  // forced biasing
  if(nullptr != biasManager) { 
    biasManager->Initialise(part,GetProcessName(),verboseLevel); 
    biasFlag = false; 
  }

  // defined ID of secondary particles
  G4int stype = GetProcessSubType();
  if(stype == fBremsstrahlung) {
    secID = _Bremsstruhlung;
    biasID = _SplitBremsstrahlung;
  } else if(stype == fPairProdByCharged) {
    secID = _PairProduction;
    mainSecondaries = 2;
  }
  baseMat = bld->GetBaseMaterialFlag();

  // initialisation of models
  numberOfModels = modelManager->NumberOfModels();
  for(G4int i=0; i<numberOfModels; ++i) {
    G4VEmModel* mod = modelManager->GetModel(i);
    if(0 == i) { currentModel = mod; }
    mod->SetMasterThread(isMaster);
    mod->SetAngularGeneratorFlag(
      theParameters->UseAngularGeneratorForIonisation());
    if(mod->HighEnergyLimit() > maxKinEnergy) {
      mod->SetHighEnergyLimit(maxKinEnergy);
    }
    mod->SetUseBaseMaterials(baseMat);
    SetEmModel(mod);
  }
  theCuts = modelManager->Initialise(particle, secondaryParticle, 
                                     1.0, verboseLevel);

  // subcut processor
  if(isIonisation) { 
    subcutProducer = lManager->SubCutProducer();
  }
  if(1 == nSCoffRegions) {
    if((*scoffRegions)[0]->GetName() == "DefaultRegionForTheWorld") {
      delete scoffRegions;
      scoffRegions = nullptr;
      nSCoffRegions = 0;
    }
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::PrepearPhysicsTable() is done "
           << " for local " << particle->GetParticleName()
           << " isIon= " << isIon;
    if(baseParticle) { 
      G4cout << "; base: " << baseParticle->GetParticleName(); 
    }
    G4cout << " chargeSqRatio= " << chargeSqRatio
           << " massRatio= " << massRatio
           << " reduceFactor= " << reduceFactor << G4endl;
    if (nSCoffRegions > 0) {
      G4cout << " SubCut secondary production is ON for regions: " << G4endl;
      for (G4int i=0; i<nSCoffRegions; ++i) {
        const G4Region* r = (*scoffRegions)[i];
        G4cout << "           " << r->GetName() << G4endl;
      }
    } else if(nullptr != subcutProducer) {
      G4cout << " SubCut secondary production is ON for all regions" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {  
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << "; local: " << particle->GetParticleName();
    if(baseParticle) { 
      G4cout << "; base: " << baseParticle->GetParticleName(); 
    }
    G4cout << " TablesAreBuilt= " << tablesAreBuilt
           << " isIon= " << isIon << "  " << this << G4endl;
  }

  if(&part == particle) {

    if(isMaster) {
      lManager->BuildPhysicsTable(particle, this);

    } else {

      const G4VEnergyLossProcess* masterProcess = 
        static_cast<const G4VEnergyLossProcess*>(GetMasterProcess());

      // copy table pointers from master thread
      SetDEDXTable(masterProcess->DEDXTable(),fRestricted);
      SetDEDXTable(masterProcess->DEDXunRestrictedTable(),fTotal);
      SetDEDXTable(masterProcess->IonisationTable(),fIsIonisation);
      SetRangeTableForLoss(masterProcess->RangeTableForLoss());
      SetCSDARangeTable(masterProcess->CSDARangeTable());
      SetSecondaryRangeTable(masterProcess->SecondaryRangeTable());
      SetInverseRangeTable(masterProcess->InverseRangeTable());
      SetLambdaTable(masterProcess->LambdaTable());
      SetTwoPeaksXS(masterProcess->TwoPeaksXS());
      isIonisation = masterProcess->IsIonisationProcess();
      baseMat = masterProcess->UseBaseMaterial();

      tablesAreBuilt = true;  
      // local initialisation of models
      G4bool printing = true;
      for(G4int i=0; i<numberOfModels; ++i) {
        G4VEmModel* mod = GetModelByIndex(i, printing);
        G4VEmModel* mod0= masterProcess->GetModelByIndex(i, printing);
        mod->SetUseBaseMaterials(baseMat);
        mod->InitialiseLocal(particle, mod0);
      }
      lManager->LocalPhysicsTables(particle, this);
    }
   
    // needs to be done only once
    safetyHelper->InitialiseHelper();
  }
  // Added tracking cut to avoid tracking artifacts
  // and identified deexcitation flag
  if(isIonisation) { 
    atomDeexcitation = lManager->AtomDeexcitation();
    if(nullptr != atomDeexcitation) { 
      if(atomDeexcitation->IsPIXEActive()) { useDeexcitation = true; } 
    }
  }

  // protection against double printout
  if(theParameters->IsPrintLocked()) { return; }

  // explicitly defined printout by particle name
  G4String num = part.GetParticleName();
  if(1 < verboseLevel || 
     (0 < verboseLevel && (num == "e-" || 
                           num == "e+"    || num == "mu+" || 
                           num == "mu-"   || num == "proton"|| 
                           num == "pi+"   || num == "pi-" || 
                           num == "kaon+" || num == "kaon-" || 
                           num == "alpha" || num == "anti_proton" || 
                           num == "GenericIon"|| num == "alpha+" )))
    { 
      StreamInfo(G4cout, part); 
    }

  /*  
  G4cout << "** G4VEnergyLossProcess::BuildPhysicsTable() for " 
         << GetProcessName() << " and " << particle->GetParticleName()
         << " isMaster: " << isMaster << " isIonisation: " 
         << isIonisation << G4endl;
  G4cout << " theDEDX: " << theDEDXTable 
         << " theRange: " << theRangeTableForLoss
         << " theInverse: " << theInverseRangeTable
         << " theLambda: " << theLambdaTable << G4endl;
  */
  //if(1 < verboseLevel || verb) {
  if(1 < verboseLevel) {
    G4cout << "### G4VEnergyLossProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName();
    if(isIonisation) { G4cout << "  isIonisation flag=1"; }
    G4cout << " baseMat=" << baseMat << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildDEDXTable(G4EmTableType tType)
{
  if(1 < verboseLevel ) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable() of type " << tType
           << " for " << GetProcessName()
           << " and particle " << particle->GetParticleName()
           << G4endl;
  }
  G4PhysicsTable* table = nullptr;
  G4double emax = maxKinEnergy;
  G4int bin = nBins;

  if(fTotal == tType) {
    emax  = maxKinEnergyCSDA;
    bin   = nBinsCSDA;
    table = theDEDXunRestrictedTable;
  } else if(fRestricted == tType) {
    table = theDEDXTable;
  } else {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable WARNING: wrong type "
           << tType << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  if(1 < verboseLevel) {
    G4cout << numOfCouples << " materials"
           << " minKinEnergy= " << minKinEnergy
           << " maxKinEnergy= " << emax
           << " nbin= " << bin
           << " EmTableType= " << tType
           << " table= " << table << "  " << this 
           << G4endl;
  }
  if(nullptr == table) { return table; }

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  G4PhysicsLogVector* aVector = nullptr;
  G4PhysicsLogVector* bVector = nullptr;

  for(size_t i=0; i<numOfCouples; ++i) {

    if(1 < verboseLevel) {
      G4cout << "G4VEnergyLossProcess::BuildDEDXVector Idx= " << i 
             << "  flagTable=  " << table->GetFlag(i) 
             << " Flag= " << bld->GetFlag(i) << G4endl;
    }
    if(bld->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple(i);
      if(nullptr != (*table)[i]) { delete (*table)[i]; }
      if(nullptr != bVector) {
        aVector = new G4PhysicsLogVector(*bVector);
      } else {
        bVector = new G4PhysicsLogVector(minKinEnergy, emax, bin, spline);
        aVector = bVector;
      }

      modelManager->FillDEDXVector(aVector, couple, tType);
      if(spline) { aVector->FillSecondDerivatives(); }

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable(): table is built for "
           << particle->GetParticleName()
           << " and process " << GetProcessName()
           << G4endl;
    if(2 < verboseLevel) G4cout << (*table) << G4endl;
  }

  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaTable(G4EmTableType tType)
{
  G4PhysicsTable* table = nullptr;

  if(fRestricted == tType) {
    table = theLambdaTable;
  } else {
    G4cout << "G4VEnergyLossProcess::BuildLambdaTable WARNING: wrong type "
           << tType << G4endl;
  }

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildLambdaTable() of type "
           << tType << " for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName()
           << " EmTableType= " << tType
           << " table= " << table
           << G4endl;
  }
  if(nullptr == table) { return table; }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  G4LossTableBuilder* bld = lManager->GetTableBuilder();

  G4PhysicsLogVector* aVector = nullptr;
  G4double scale = G4Log(maxKinEnergy/minKinEnergy);

  for(size_t i=0; i<numOfCouples; ++i) {

    if (bld->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = 
        theCoupleTable->GetMaterialCutsCouple(i);
      delete (*table)[i];

      G4bool startNull = true;
      G4double emin = 
        MinPrimaryEnergy(particle,couple->GetMaterial(),(*theCuts)[i]);
      if(minKinEnergy > emin) { 
        emin = minKinEnergy; 
        startNull = false;
      }

      G4double emax = maxKinEnergy;
      if(emax <= emin) { emax = 2*emin; }
      G4int bin = G4lrint(nBins*G4Log(emax/emin)/scale);
      bin = std::max(bin, 3);
      aVector = new G4PhysicsLogVector(emin, emax, bin, spline);

      modelManager->FillLambdaVector(aVector, couple, startNull, tType);
      if(spline) { aVector->FillSecondDerivatives(); }

      // Insert vector for this material into the table
      G4PhysicsTableHelper::SetPhysicsVector(table, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built for "
           << particle->GetParticleName()
           << G4endl;
  }

  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::StreamInfo(std::ostream& out,
                const G4ParticleDefinition& part, G4bool rst) const
{
  G4String indent = (rst ? "  " : "");
  out << std::setprecision(6);
  out << G4endl << indent << GetProcessName()  << ": ";
  if (!rst) out << " for " << part.GetParticleName();
  out << "  XStype:" << fXSType 
      << "  SubType=" << GetProcessSubType() << G4endl
      << "      dE/dx and range tables from "
      << G4BestUnit(minKinEnergy,"Energy")
      << " to " << G4BestUnit(maxKinEnergy,"Energy")
      << " in " << nBins << " bins" << G4endl
      << "      Lambda tables from threshold to "
      << G4BestUnit(maxKinEnergy,"Energy")
      << ", " << theParameters->NumberOfBinsPerDecade() 
      << " bins/decade, spline: " << spline
      << G4endl;
  if(nullptr != theRangeTableForLoss && isIonisation) {
    out << "      StepFunction=(" << dRoverRange << ", "
        << finalRange/mm << " mm)"
        << ", integ: " << fXSType
        << ", fluct: " << lossFluctuationFlag
        << ", linLossLim= " << linLossLimit
        << G4endl;
  }
  StreamProcessInfo(out);
  modelManager->DumpModelList(out, verboseLevel);
  if(nullptr != theCSDARangeTable && isIonisation) {
    out << "      CSDA range table up"
        << " to " << G4BestUnit(maxKinEnergyCSDA,"Energy")
        << " in " << nBinsCSDA << " bins" << G4endl;
  }
  if(nSCoffRegions>0 && isIonisation) {
    out << "      Subcutoff sampling in " << nSCoffRegions 
        << " regions" << G4endl;
  }
  if(2 < verboseLevel) {
    out << "      DEDXTable address= " << theDEDXTable << G4endl; 
    if(nullptr != theDEDXTable && isIonisation) 
      out << (*theDEDXTable) << G4endl;
    out << "non restricted DEDXTable address= " 
        << theDEDXunRestrictedTable << G4endl;
    if(nullptr != theDEDXunRestrictedTable && isIonisation) {
      out << (*theDEDXunRestrictedTable) << G4endl;
    }
    out << "      CSDARangeTable address= " << theCSDARangeTable << G4endl;
    if(nullptr != theCSDARangeTable && isIonisation) {
      out << (*theCSDARangeTable) << G4endl;
    }
    out << "      RangeTableForLoss address= " << theRangeTableForLoss 
        << G4endl;
    if(nullptr != theRangeTableForLoss && isIonisation) {
      out << (*theRangeTableForLoss) << G4endl;
    }
    out << "      InverseRangeTable address= " << theInverseRangeTable 
        << G4endl;
    if(nullptr != theInverseRangeTable && isIonisation) {
      out << (*theInverseRangeTable) << G4endl;
    }
    out << "      LambdaTable address= " << theLambdaTable << G4endl;
    if(nullptr != theLambdaTable) {
      out << (*theLambdaTable) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ActivateSubCutoff(const G4Region* r)
{
  if(nullptr == scoffRegions) {
    scoffRegions = new std::vector<const G4Region*>;
  }
  // the region is in the list
  if(!scoffRegions->empty()) {
    for (auto & reg : *scoffRegions) {
      if (reg == r) { return; }
    }
  }
  // new region 
  scoffRegions->push_back(r);
  ++nSCoffRegions;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::IsRegionForCubcutProcessor(const G4Track& aTrack)
{
  if(0 == nSCoffRegions) { return true; }
  const G4Region* r = aTrack.GetVolume()->GetLogicalVolume()->GetRegion();
  for(auto & reg : *scoffRegions) {
    if(r == reg) { return true; }
  }
  return false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::StartTracking(G4Track* track)
{
  /*
  G4cout << "G4VEnergyLossProcess::StartTracking: " 
         << track->GetDefinition()->GetParticleName() 
         << " e(MeV)= " << track->GetKineticEnergy();
  if(particle) G4cout << "  " << particle->GetParticleName();
  if(baseParticle) G4cout << " basePart: " << baseParticle->GetParticleName();
  G4cout << "  " << GetProcessName();
  if(isIon) G4cout << " isIon:  Q=" << track->GetDefinition()->GetPDGCharge() 
  << " Qdyn=" << track->GetDynamicParticle()->GetCharge(); 
  G4cout << G4endl;
  */
  // reset parameters for the new track
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX;
  currentCouple = nullptr;

  // reset ion
  if(isIon) {
    const G4double newmass = track->GetDefinition()->GetPDGMass();
    if(nullptr != baseParticle) {
      massRatio    = baseParticle->GetPDGMass()/newmass;
      logMassRatio = G4Log(massRatio);
    } else if(nullptr != theGenericIon) {
      massRatio    = CLHEP::proton_mass_c2/newmass;
      logMassRatio = G4Log(massRatio);
    } else {
      massRatio    = 1.0;
      logMassRatio = 0.0;
    }
  }  
  // forced biasing only for primary particles
  if(biasManager) {
    if(0 == track->GetParentID()) {
      biasFlag = true; 
      biasManager->ResetForcedInteraction(); 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::AlongStepGetPhysicalInteractionLength(
                             const G4Track&,G4double,G4double,G4double&,
                             G4GPILSelection* selection)
{
  G4double x = DBL_MAX;
  *selection = aGPILSelection;
  if(isIonisation && currentModel->IsActive(preStepScaledEnergy)) {
    GetScaledRangeForScaledEnergy(preStepScaledEnergy, preStepLogScaledEnergy);
    const G4double finR = (rndmStepFlag) ? std::min(finalRange,
      currentCouple->GetProductionCuts()->GetProductionCut(1)) : finalRange;
    x = (fRange > finR) ? 
      fRange*dRoverRange + finR*(1.0-dRoverRange)*(2.0-finR/fRange) : fRange; 
    // if(particle->GetPDGMass() > 0.9*GeV)
    /*
    G4cout << GetProcessName() << ": e= " << preStepKinEnergy
           <<" range= "<<fRange << " idx= " << basedCoupleIndex
           << " finR= " << finR << " limit= " << x <<
           << "\n" << "massRatio= " << massRatio << " Q^2= " << chargeSqRatio 
           << " dRoverRange= " << dRoverRange 
           << " finalRange= " << finalRange << G4endl;
    */
  }
  //G4cout<<"AlongStepGPIL: " << GetProcessName()<<": e= "<<preStepKinEnergy 
  //<<" stepLimit= "<<x<<G4endl;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  G4double x = DBL_MAX;

  // initialisation of material, mass, charge, model 
  // at the beginning of the step
  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy       = track.GetKineticEnergy();
  preStepLogKinEnergy    = track.GetDynamicParticle()->GetLogKineticEnergy();
  preStepScaledEnergy    = preStepKinEnergy*massRatio;
  preStepLogScaledEnergy = preStepLogKinEnergy + logMassRatio;
  SelectModel(preStepScaledEnergy);

  if(!currentModel->IsActive(preStepScaledEnergy)) { 
    theNumberOfInteractionLengthLeft = -1.0;
    currentInteractionLength = DBL_MAX;
    return x; 
  }

  // change effective charge of a charged particle on fly
  if(isIon) {
    const G4double q2 = currentModel->ChargeSquareRatio(track);
    if(q2 != chargeSqRatio) { 
      fFactor *= q2/chargeSqRatio;
      reduceFactor = 1.0/(fFactor*massRatio);
      chargeSqRatio = q2;
      // G4cout << "PostStepGPIL: Q^2=" << chargeSqRatio << " reducedFactor=" << reduceFactor << G4endl;
    }
  }

  // forced biasing only for primary particles
  if(biasManager) {
    if(0 == track.GetParentID() && biasFlag && 
       biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
      return biasManager->GetStepLimit(currentCoupleIndex, previousStepSize);
    }
  }

  // compute mean free path
  ComputeLambdaForScaledEnergy(preStepScaledEnergy, preStepLogScaledEnergy);

  // zero cross section
  if(preStepLambda <= 0.0) { 
    theNumberOfInteractionLengthLeft = -1.0;
    currentInteractionLength = DBL_MAX;
  } else {

    // non-zero cross section
    if (theNumberOfInteractionLengthLeft < 0.0) {

      // beggining of tracking (or just after DoIt of this process)
      theNumberOfInteractionLengthLeft = -G4Log( G4UniformRand() );
      theInitialNumberOfInteractionLength = theNumberOfInteractionLengthLeft; 

    } else if(currentInteractionLength < DBL_MAX) {

      // subtract NumberOfInteractionLengthLeft using previous step
      theNumberOfInteractionLengthLeft -= 
        previousStepSize/currentInteractionLength;

      theNumberOfInteractionLengthLeft = 
        std::max(theNumberOfInteractionLengthLeft, 0.0);
    }

    // new mean free path and step limit
    currentInteractionLength = 1.0/preStepLambda;
    x = theNumberOfInteractionLengthLeft * currentInteractionLength;
  }
#ifdef G4VERBOSE
  if (verboseLevel>2){
    //  if(particle->GetPDGMass() > 0.9*GeV){
    G4cout << "G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" << G4endl; 
    G4cout << " for " << track.GetDefinition()->GetParticleName() 
           << " in Material  " <<  currentMaterial->GetName()
           << " Ekin(MeV)= " << preStepKinEnergy/MeV 
           << "  " << track.GetMaterial()->GetName()
           <<G4endl;
    G4cout << "MeanFreePath = " << currentInteractionLength/cm << "[cm]" 
           << "InteractionLength= " << x/cm <<"[cm] " <<G4endl;
  }
#endif
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void
G4VEnergyLossProcess::ComputeLambdaForScaledEnergy(G4double e, G4double loge)
{
  // cross section increased with energy
  if(fXSType == fEmIncreasing) {
    if(e/lambdaFactor < mfpKinEnergy) {
      mfpKinEnergy = e;
      preStepLambda = GetLambdaForScaledEnergy(e, loge); 
    }

  // cross section has two peaks
  } else if(fXSType == fEmTwoPeaks) {
    G4TwoPeaksXS* xs = (*fXSpeaks)[basedCoupleIndex];
    const G4double e1peak = xs->e1peak;

    // below the 1st peak
    if(e <= e1peak) {
      if(e/lambdaFactor < mfpKinEnergy) {
        mfpKinEnergy = e;
        preStepLambda = GetLambdaForScaledEnergy(e, loge); 
      }
      return;
    }
    const G4double e1deep = xs->e1deep;
    // above the 1st peak, below the deep
    if(e <= e1deep) {
      if(mfpKinEnergy >= e1deep || e <= mfpKinEnergy) { 
        const G4double e1 = std::max(e1peak, e*lambdaFactor);
        preStepLambda = GetLambdaForScaledEnergy(e1); 
        mfpKinEnergy = e1;
      }
      return;
    }
    const G4double e2peak = xs->e2peak;
    // above the deep, below 2nd peak
    if(e <= e2peak) {
      if(e/lambdaFactor < mfpKinEnergy) {
        mfpKinEnergy = e;
        preStepLambda = GetLambdaForScaledEnergy(e, loge); 
      }
      return;
    }
    const G4double e2deep = xs->e2deep;
    // above the 2nd peak, below the deep
    if(e <= e2deep) {
      if(mfpKinEnergy >= e2deep || e <= mfpKinEnergy) { 
        const G4double e1 = std::max(e2peak, e*lambdaFactor);
        preStepLambda = GetLambdaForScaledEnergy(e1); 
        mfpKinEnergy = e1;
      }
      return;
    }
    // above the deep, below 3d peak
    if(e/lambdaFactor < mfpKinEnergy) {
      mfpKinEnergy = e;
      preStepLambda = GetLambdaForScaledEnergy(e, loge); 
    }

    // integral method is not used
  } else {
    preStepLambda = GetLambdaForScaledEnergy(e, loge); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::AlongStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
  fParticleChange.InitializeForAlongStep(track);
  // The process has range table - calculate energy loss
  if(!isIonisation || !currentModel->IsActive(preStepScaledEnergy)) {
    return &fParticleChange;
  }

  // Get the actual (true) Step length
  G4double length = step.GetStepLength();
  if(length <= 0.0) { return &fParticleChange; }
  G4double eloss  = 0.0;
 
  /*
  if(-1 < verboseLevel) {
    const G4ParticleDefinition* d = track.GetParticleDefinition();
    G4cout << "AlongStepDoIt for "
           << GetProcessName() << " and particle "
           << d->GetParticleName()
           << "  eScaled(MeV)= " << preStepScaledEnergy/MeV
           << "  range(mm)= " << fRange/mm
           << "  s(mm)= " << length/mm
           << "  rf= " << reduceFactor
           << "  q^2= " << chargeSqRatio
           << " md= " << d->GetPDGMass()
           << "  status= " << track.GetTrackStatus()
           << "  " << track.GetMaterial()->GetName()
           << G4endl;
  }
  */

  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();
  if(weightFlag) {
    weight /= biasFactor;
    fParticleChange.ProposeWeight(weight);
  }

  // stopping
  if (length >= fRange || preStepKinEnergy <= lowestKinEnergy) {
    eloss = preStepKinEnergy;
    if (useDeexcitation) {
      atomDeexcitation->AlongStepDeexcitation(scTracks, step, 
                                              eloss, currentCoupleIndex);
      if(scTracks.size() > 0) { FillSecondariesAlongStep(weight); }
      eloss = std::max(eloss, 0.0);
    }
    fParticleChange.SetProposedKineticEnergy(0.0);
    fParticleChange.ProposeLocalEnergyDeposit(eloss);
    return &fParticleChange;
  }
  //G4cout << theDEDXTable << "  idx= " << basedCoupleIndex 
  // << "  " << GetProcessName() << "  "<< currentMaterial->GetName()<<G4endl;
  //if(particle->GetParticleName() == "e-")G4cout << (*theDEDXTable) <<G4endl;
  // Short step
  eloss = length*GetDEDXForScaledEnergy(preStepScaledEnergy,
                                        preStepLogScaledEnergy);

  //G4cout << "Short STEP: eloss= " << eloss << G4endl;

  // Long step
  if(eloss > preStepKinEnergy*linLossLimit) {

    G4double x = (fRange - length)/reduceFactor;
    //G4cout << "x= " << x << "  " << theInverseRangeTable << G4endl;
    eloss = preStepKinEnergy - ScaledKinEnergyForLoss(x)/massRatio;
   
    /*
    if(-1 < verboseLevel) 
      G4cout << "Long STEP: rPre(mm)= " 
             << GetScaledRangeForScaledEnergy(preStepScaledEnergy)/mm
             << " rPost(mm)= " << x/mm
             << " ePre(MeV)= " << preStepScaledEnergy/MeV
             << " eloss(MeV)= " << eloss/MeV
             << " eloss0(MeV)= "
             << GetDEDXForScaledEnergy(preStepScaledEnergy)*length/MeV
             << " lim(MeV)= " << preStepKinEnergy*linLossLimit/MeV
             << G4endl;
    */
  }

  /*
  G4double eloss0 = eloss;
  if(-1 < verboseLevel ) {
    G4cout << "Before fluct: eloss(MeV)= " << eloss/MeV
           << " e-eloss= " << preStepKinEnergy-eloss
           << " step(mm)= " << length/mm
           << " range(mm)= " << fRange/mm
           << " fluct= " << lossFluctuationFlag
           << G4endl;
  }
  */

  const G4double cut = (*theCuts)[currentCoupleIndex];
  G4double esec = 0.0;

  // Corrections, which cannot be tabulated
  if(isIon) {
    currentModel->CorrectionsAlongStep(currentCouple, dynParticle, 
                                       length, eloss);
    eloss = std::max(eloss, 0.0);
  }

  // Sample fluctuations if not full energy loss
  if(eloss >= preStepKinEnergy) {
    eloss = preStepKinEnergy;

  } else if (lossFluctuationFlag) {
    const G4double tmax = currentModel->MaxSecondaryKinEnergy(dynParticle);
    const G4double tcut = std::min(cut, tmax);
    G4VEmFluctuationModel* fluc = currentModel->GetModelOfFluctuations();
    eloss = fluc->SampleFluctuations(currentCouple,dynParticle,
                                     tcut, tmax, length, eloss);
    /*
    if(-1 < verboseLevel) 
      G4cout << "After fluct: eloss(MeV)= " << eloss/MeV
             << " fluc= " << (eloss-eloss0)/MeV
             << " ChargeSqRatio= " << chargeSqRatio
             << " massRatio= " << massRatio
             << " tmax= " << tmax
             << G4endl;
    */
  }

  // deexcitation
  if (useDeexcitation) {
    G4double esecfluo = preStepKinEnergy;
    G4double de = esecfluo;
    //G4double eloss0 = eloss;
    /*
    G4cout << "### 1: E(keV)= " << preStepKinEnergy/keV
           << " Efluomax(keV)= " << de/keV
           << " Eloss(keV)= " << eloss/keV << G4endl; 
    */
    atomDeexcitation->AlongStepDeexcitation(scTracks, step, 
                                            de, currentCoupleIndex);

    // sum of de-excitation energies
    esecfluo -= de;

    // subtracted from energy loss
    if(eloss >= esecfluo) {
      esec  += esecfluo;
      eloss -= esecfluo;
    } else {
      esec += esecfluo;
      eloss = 0.0; 
    } 
    /*    
    if(esecfluo > 0.0) {
      G4cout << "### 2: E(keV)= " << preStepKinEnergy/keV
             << " Esec(keV)= " << esec/keV
             << " Esecf(kV)= " << esecfluo/keV
             << " Eloss0(kV)= " << eloss0/keV
             << " Eloss(keV)= " << eloss/keV 
             << G4endl; 
    } 
    */   
  }
  if(nullptr != subcutProducer && IsRegionForCubcutProcessor(track)) {
    subcutProducer->SampleSecondaries(step, scTracks, eloss, cut);
  }
  // secondaries from atomic de-excitation and subcut
  if(!scTracks.empty()) { FillSecondariesAlongStep(weight); }

  // Energy balance
  G4double finalT = preStepKinEnergy - eloss - esec;
  if (finalT <= lowestKinEnergy) {
    eloss += finalT;
    finalT = 0.0;
  } else if(isIon) {
    fParticleChange.SetProposedCharge(
      currentModel->GetParticleCharge(track.GetParticleDefinition(),
                                      currentMaterial,finalT));
  }
  eloss = std::max(eloss, 0.0);

  fParticleChange.SetProposedKineticEnergy(finalT);
  fParticleChange.ProposeLocalEnergyDeposit(eloss);
  /*
  if(-1 < verboseLevel) {
    G4double del = finalT + eloss + esec - preStepKinEnergy;
    G4cout << "Final value eloss(MeV)= " << eloss/MeV
           << " preStepKinEnergy= " << preStepKinEnergy
           << " postStepKinEnergy= " << finalT
           << " de(keV)= " << del/keV
           << " lossFlag= " << lossFluctuationFlag
           << "  status= " << track.GetTrackStatus()
           << G4endl;
  }
  */
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::FillSecondariesAlongStep(G4double wt)
{
  const G4int n0 = scTracks.size();
  G4double weight = wt;
  // weight may be changed by biasing manager
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion(currentCoupleIndex)) {
      weight *=
        biasManager->ApplySecondaryBiasing(scTracks, currentCoupleIndex);
    }
  } 

  // fill secondaries
  const G4int n = scTracks.size();
  fParticleChange.SetNumberOfSecondaries(n);

  for(G4int i=0; i<n; ++i) {
    G4Track* t = scTracks[i];
    if(nullptr != t) {
      t->SetWeight(weight); 
      pParticleChange->AddSecondary(t);
      if(i >= n0) { t->SetCreatorModelID(biasID); }
      //G4cout << "Secondary(along step) has weight " << t->GetWeight() 
      //<< ", kenergy " << t->GetKineticEnergy()/MeV << " MeV" <<G4endl;
    }
  }
  scTracks.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::PostStepDoIt(const G4Track& track,
                                                      const G4Step& step)
{
  // In all cases clear number of interaction lengths
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX;

  fParticleChange.InitializeForPostStep(track);
  const G4double finalT = track.GetKineticEnergy();

  const G4double postStepScaledEnergy = finalT*massRatio;
  SelectModel(postStepScaledEnergy);

  if(!currentModel->IsActive(postStepScaledEnergy)) { 
    return &fParticleChange; 
  }
  /*
  if(-1 < verboseLevel) {
    G4cout << GetProcessName()
           << "::PostStepDoIt: E(MeV)= " << finalT/MeV
           << G4endl;
  }
  */

  // forced process - should happen only once per track
  if(biasFlag) {
    if(biasManager->ForcedInteractionRegion(currentCoupleIndex)) {
      biasFlag = false;
    }
  }

  const G4DynamicParticle* dp = track.GetDynamicParticle();

  // Integral approach
  if (fXSType != fEmNoIntegral) {
    const G4double logFinalT = dp->GetLogKineticEnergy();
    G4double lx = GetLambdaForScaledEnergy(postStepScaledEnergy,
                                           logFinalT + logMassRatio);
    lx = std::max(lx, 0.0);

    // cache cross section useful for the false interaction
    const G4double lg = preStepLambda;
    if(postStepScaledEnergy < mfpKinEnergy) {
      mfpKinEnergy = postStepScaledEnergy;
      preStepLambda = lx;
    }
    /*
    if(preStepLambda<lx && 1 < verboseLevel) {
      G4cout << "WARNING: for " << particle->GetParticleName()
             << " and " << GetProcessName()
             << " E(MeV)= " << finalT/MeV
             << " preLambda= " << preStepLambda 
             << " < " << lx << " (postLambda) "
             << G4endl;
    }
    */
    // if both lg and lx are zero then no interaction
    if(lg*G4UniformRand() >= lx) {
      return &fParticleChange;
    }
  }
  // define new weight for primary and secondaries
  G4double weight = fParticleChange.GetParentWeight();
  if(weightFlag) {
    weight /= biasFactor;
    fParticleChange.ProposeWeight(weight);
  }

  const G4double tcut = (*theCuts)[currentCoupleIndex];

  // sample secondaries
  secParticles.clear();
  //G4cout<< "@@@ Eprimary= "<<dynParticle->GetKineticEnergy()/MeV
  //        << " cut= " << tcut/MeV << G4endl;
  currentModel->SampleSecondaries(&secParticles, currentCouple, dp, tcut);

  const G4int num0 = secParticles.size();

  // bremsstrahlung splitting or Russian roulette  
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion(currentCoupleIndex)) {
      G4double eloss = 0.0;
      weight *= biasManager->ApplySecondaryBiasing(
                                      secParticles,
                                      track, currentModel, 
                                      &fParticleChange, eloss,
                                      currentCoupleIndex, tcut, 
                                      step.GetPostStepPoint()->GetSafety());
      if(eloss > 0.0) {
        eloss += fParticleChange.GetLocalEnergyDeposit();
        fParticleChange.ProposeLocalEnergyDeposit(eloss);
      }
    }
  }

  // save secondaries
  const G4int num = secParticles.size();
  if(num > 0) {

    fParticleChange.SetNumberOfSecondaries(num);
    G4double time = track.GetGlobalTime();

    G4int n1(0), n2(0);
    if(num0 > mainSecondaries) { 
      currentModel->FillNumberOfSecondaries(n1, n2);
    }

    for (G4int i=0; i<num; ++i) {
      if(nullptr != secParticles[i]) {
        G4Track* t = new G4Track(secParticles[i], time, track.GetPosition());
        t->SetTouchableHandle(track.GetTouchableHandle());
        if (biasManager) {
          t->SetWeight(weight * biasManager->GetWeight(i));
        } else {
          t->SetWeight(weight);
        }
        if(i < num0) {
          t->SetCreatorModelID(secID);
        } else if(i < num0 + n1) {
          t->SetCreatorModelID(tripletID);
        } else {
          t->SetCreatorModelID(biasID);
        }

        //G4cout << "Secondary(post step) has weight " << t->GetWeight() 
        //       << ", kenergy " << t->GetKineticEnergy()/MeV << " MeV" 
        //       << " time= " << time/ns << " ns " << G4endl;
        pParticleChange->AddSecondary(t);
      }
    }
  }

  if(0.0 == fParticleChange.GetProposedKineticEnergy() &&
     fAlive == fParticleChange.GetTrackStatus()) {
    if(particle->GetProcessManager()->GetAtRestProcessVector()->size() > 0)
         { fParticleChange.ProposeTrackStatus(fStopButAlive); }
    else { fParticleChange.ProposeTrackStatus(fStopAndKill); }
  }

  /*
  if(-1 < verboseLevel) {
    G4cout << "::PostStepDoIt: Sample secondary; Efin= " 
    << fParticleChange.GetProposedKineticEnergy()/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit()
           << ", " <<  currentModel->HighEnergyLimit() << ")"
           << "  preStepLambda= " << preStepLambda
           << "  dir= " << track.GetMomentumDirection()
           << "  status= " << track.GetTrackStatus()
           << G4endl;
  }
  */
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEnergyLossProcess::StorePhysicsTable(
       const G4ParticleDefinition* part, const G4String& directory, 
       G4bool ascii)
{
  G4bool res = true;
  //G4cout << "G4VEnergyLossProcess::StorePhysicsTable: " << part->GetParticleName()
  //         << "  " << directory << "  " << ascii << G4endl;
  if (!isMaster || baseParticle || part != particle ) return res;

  if(!StoreTable(part,theDEDXTable,ascii,directory,"DEDX")) 
    {res = false;}

  if(!StoreTable(part,theDEDXunRestrictedTable,ascii,directory,"DEDXnr")) 
    {res = false;}

  if(!StoreTable(part,theIonisationTable,ascii,directory,"Ionisation")) 
    {res = false;}

  if(isIonisation &&
     !StoreTable(part,theCSDARangeTable,ascii,directory,"CSDARange")) 
    {res = false;}

  if(isIonisation &&
     !StoreTable(part,theRangeTableForLoss,ascii,directory,"Range")) 
    {res = false;}
  
  if(isIonisation &&
     !StoreTable(part,theInverseRangeTable,ascii,directory,"InverseRange")) 
    {res = false;}
  
  if(!StoreTable(part,theLambdaTable,ascii,directory,"Lambda")) 
    {res = false;}

  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool 
G4VEnergyLossProcess::RetrievePhysicsTable(const G4ParticleDefinition* part, 
                                           const G4String& directory,
                                           G4bool ascii)
{
  G4bool res = true;
  if (!isMaster) return res;
  const G4String& particleName = part->GetParticleName();

  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::RetrievePhysicsTable() for "
           << particleName << " and process " << GetProcessName()
           << "; tables_are_built= " << tablesAreBuilt
           << G4endl;
  }
  if(particle == part) {

    if(nullptr == baseParticle) {

      G4bool fpi = true;
      if(!RetrieveTable(part,theDEDXTable,ascii,directory,"DEDX",fpi)) 
        { fpi = false; }

      // ionisation table keeps individual dEdx and not sum of sub-processes
      if(!RetrieveTable(part,theDEDXTable,ascii,directory,"Ionisation",false)) 
        { fpi = false; }

      if(!RetrieveTable(part,theRangeTableForLoss,ascii,directory,"Range",fpi)) 
        { res = false; }

      if(!RetrieveTable(part,theDEDXunRestrictedTable,ascii,directory,
                        "DEDXnr",false)) 
        { res = false; }

      if(!RetrieveTable(part,theCSDARangeTable,ascii,directory,
                        "CSDARange",false)) 
        { res = false; }

      if(!RetrieveTable(part,theInverseRangeTable,ascii,directory,
                        "InverseRange",fpi)) 
        { res = false; }

      if(!RetrieveTable(part,theLambdaTable,ascii,directory,"Lambda",true)) 
        { res = false; }
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool G4VEnergyLossProcess::StoreTable(const G4ParticleDefinition* part, 
                                        G4PhysicsTable* aTable, G4bool ascii,
                                        const G4String& directory,
                                        const G4String& tname)
{
  G4bool res = true;
  if (nullptr != aTable) {
    const G4String& name = GetPhysicsTableFileName(part, directory, tname, ascii);
    if ( aTable->StorePhysicsTable(name,ascii) ) {
      if (0 < verboseLevel) G4cout << "Stored: " << name << G4endl;
    } else {
      res = false;
      G4cout << "Fail to store: " << name << G4endl;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool 
G4VEnergyLossProcess::RetrieveTable(const G4ParticleDefinition* part, 
                                    G4PhysicsTable* aTable, 
                                    G4bool ascii,
                                    const G4String& directory,
                                    const G4String& tname,
                                    G4bool mandatory)
{
  G4bool isRetrieved = false;
  G4String filename = GetPhysicsTableFileName(part,directory,tname,ascii);
  if(nullptr != aTable) {
    if(aTable->ExistPhysicsTable(filename)) {
      if(G4PhysicsTableHelper::RetrievePhysicsTable(aTable,filename,ascii,spline)) {
        isRetrieved = true;
        if(spline) {
          for(auto & v : *aTable) { 
            if(nullptr != v) { v->FillSecondDerivatives(); } 
          }
        }
        if (0 < verboseLevel) {
          G4cout << tname << " table for " << part->GetParticleName() 
                 << " is Retrieved from <" << filename << ">"
                 << G4endl;
        }
      }
    }
  }
  if(mandatory && !isRetrieved) {
    if(0 < verboseLevel) {
      G4cout << tname << " table for " << part->GetParticleName() 
             << " from file <"
             << filename << "> is not Retrieved"
             << G4endl;
    }
    return false;
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::GetDEDXDispersion(
                                  const G4MaterialCutsCouple *couple,
                                  const G4DynamicParticle* dp,
                                        G4double length)
{
  DefineMaterial(couple);
  G4double ekin = dp->GetKineticEnergy();
  SelectModel(ekin*massRatio);
  G4double tmax = currentModel->MaxSecondaryKinEnergy(dp);
  G4double tcut = std::min(tmax,(*theCuts)[currentCoupleIndex]);
  G4double d = 0.0;
  G4VEmFluctuationModel* fm = currentModel->GetModelOfFluctuations();
  if(nullptr != fm) { d = fm->Dispersion(currentMaterial,dp,tcut,tmax,length); }
  return d;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double
G4VEnergyLossProcess::CrossSectionPerVolume(G4double kineticEnergy,
                                            const G4MaterialCutsCouple* couple,
                                            G4double logKineticEnergy)
{
  // Cross section per volume is calculated
  DefineMaterial(couple);
  G4double cross = 0.0;
  if (nullptr != theLambdaTable) {
    cross = GetLambdaForScaledEnergy(kineticEnergy * massRatio,
                                     logKineticEnergy + logMassRatio);
  } else {
    SelectModel(kineticEnergy*massRatio);
    cross = (!baseMat) ? biasFactor : biasFactor*(*theDensityFactor)[currentCoupleIndex];
    cross *= (currentModel->CrossSectionPerVolume(currentMaterial, particle, kineticEnergy,
                                                  (*theCuts)[currentCoupleIndex]));
  }
  return std::max(cross, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::MeanFreePath(const G4Track& track)
{
  DefineMaterial(track.GetMaterialCutsCouple());
  const G4double kinEnergy    = track.GetKineticEnergy();
  const G4double logKinEnergy = track.GetDynamicParticle()->GetLogKineticEnergy();
  const G4double cs = GetLambdaForScaledEnergy(kinEnergy * massRatio, 
                                               logKinEnergy + logMassRatio);
  return (0.0 < cs) ? 1.0/cs : DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::ContinuousStepLimit(const G4Track& track, 
                                                   G4double x, G4double y, 
                                                   G4double& z)
{
  return AlongStepGetPhysicalInteractionLength(track, x, y, z, &aGPILSelection);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::GetMeanFreePath(
                             const G4Track& track,
                             G4double,
                             G4ForceCondition* condition)

{
  *condition = NotForced;
  return MeanFreePath(track);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEnergyLossProcess::GetContinuousStepLimit(
                const G4Track&,
                G4double, G4double, G4double&)
{
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* 
G4VEnergyLossProcess::LambdaPhysicsVector(const G4MaterialCutsCouple* couple, 
                                          G4double)
{
  DefineMaterial(couple);
  G4PhysicsVector* v = (*theLambdaTable)[basedCoupleIndex];
  return new G4PhysicsVector(*v);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::SetDEDXTable(G4PhysicsTable* p, G4EmTableType tType)
{
  if(fTotal == tType) {
    theDEDXunRestrictedTable = p;

  } else if(fRestricted == tType) {
    /*
      G4cout<< "G4VEnergyLossProcess::SetDEDXTable "
            << particle->GetParticleName()
            << " oldTable " << theDEDXTable << " newTable " << p 
            << " ion " << theIonisationTable 
            << " IsMaster " << isMaster 
            << " " << GetProcessName() << G4endl;
      G4cout << (*p) << G4endl;
    */
    theDEDXTable = p;
  } else if(fIsIonisation == tType) {
    /*
      G4cout<< "G4VEnergyLossProcess::SetIonisationTable "
            << particle->GetParticleName()
            << " oldTable " << theDEDXTable << " newTable " << p 
            << " ion " << theIonisationTable 
            << " IsMaster " << isMaster 
            << " " << GetProcessName() << G4endl;
    */
    theIonisationTable = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetCSDARangeTable(G4PhysicsTable* p)
{
  theCSDARangeTable = p; 
  if(1 < verboseLevel) {
    G4cout << "### Set CSDA Range table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetRangeTableForLoss(G4PhysicsTable* p)
{
  theRangeTableForLoss = p;
  if(1 < verboseLevel) {
    G4cout << "### Set Range table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetSecondaryRangeTable(G4PhysicsTable* p)
{
  theSecondaryRangeTable = p;
  if(1 < verboseLevel) {
    G4cout << "### Set SecondaryRange table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetInverseRangeTable(G4PhysicsTable* p)
{
  theInverseRangeTable = p;
  if(1 < verboseLevel) {
    G4cout << "### Set InverseRange table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaTable(G4PhysicsTable* p)
{
  if(1 < verboseLevel) {
    G4cout << "### Set Lambda table " << p 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
    //G4cout << *p << G4endl;
  }
  theLambdaTable = p; 
  tablesAreBuilt = true;

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  theDensityFactor = bld->GetDensityFactors();
  theDensityIdx = bld->GetCoupleIndexes();

  if(isMaster && nullptr == baseParticle && 
     nullptr != theLambdaTable && fEmTwoPeaks == fXSType) {

    size_t n = theLambdaTable->length();

    G4double e, ss, xs, ee, e1peak, xs1peak, e1deep, e2peak, e2deep, xs2peak;

    // first loop on existing vectors
    for (size_t i=0; i<n; ++i) {
      const G4PhysicsVector* pv = (*theLambdaTable)[i];
      ee = xs = xs1peak = xs2peak = 0.0;
      e1peak = e1deep = e2peak = e2deep = DBL_MAX;
      if(nullptr != pv) {
        size_t nb = pv->GetVectorLength();
        for (size_t j=0; j<nb; ++j) {
          e = pv->Energy(j);
          ss = (*pv)(j);
          // find out 1st peak
          if(e1peak == DBL_MAX) {
            if(ss >= xs) {
              xs = ss;
              ee = e;
              continue;
            } else {
              e1peak = ee;
              xs1peak = xs;
            }
          }
          // find out the deep
          if(e1deep == DBL_MAX) {
            if(ss <= xs) {
              xs = ss;
              ee = e;
              continue;
            } else {
              e1deep = ee;
            }
          }
          // find out 2nd peak
          if(e2peak == DBL_MAX) {
            if(ss >= xs) {
              xs = ss;
              ee = e;
              continue;
            } else {
              e2peak = ee;
              xs2peak = xs;
            }
          }
          if(e2deep == DBL_MAX) {
            if(ss <= xs) {
              xs = ss;
              ee = e;
              continue;
            } else {
              e2deep = ee;
              break;
            }
          }
        }
      }
      G4TwoPeaksXS* x = (*fXSpeaks)[i];
      if(nullptr == x) { 
        x = new G4TwoPeaksXS(); 
        (*fXSpeaks)[i] = x;
      }
      x->e1peak = e1peak;
      x->e1deep = e1deep;
      x->e2peak = e2peak;
      x->e2deep = e2deep;
       
      if(1 < verboseLevel) {
        G4cout << "For " << particle->GetParticleName() 
               << " index= " << i << " data:\n" << " E1peak=" << e1peak 
               << " xs1= " << xs1peak << " E1deep=" << e1deep
               << " E2peak=" << e2peak << " xs2=" << xs2peak 
               << " E2deep=" << e2deep << G4endl;
      }
    }
    // second loop using base materials
    for (size_t i=0; i<n; ++i) {
      const G4PhysicsVector* pv = (*theLambdaTable)[i];
      if (nullptr == pv) {
        G4int j = (*theDensityIdx)[i];
        G4TwoPeaksXS* x = (*fXSpeaks)[i];
        G4TwoPeaksXS* y = (*fXSpeaks)[j];
        if(nullptr == x) { 
          x = new G4TwoPeaksXS(); 
          (*fXSpeaks)[i] = x;
        }
        x->e1peak = y->e1peak;
        x->e1deep = y->e1deep;
        x->e2peak = y->e2peak;
        x->e2deep = y->e2deep;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetTwoPeaksXS(std::vector<G4TwoPeaksXS*>* ptr)
{
  fXSpeaks = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4VEnergyLossProcess::GetCurrentElement() const
{
  return (nullptr != currentModel) ? currentModel->GetCurrentElement() : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetCrossSectionBiasingFactor(G4double f, 
                                                        G4bool flag)
{
  if(f > 0.0) { 
    biasFactor = f; 
    weightFlag = flag;
    if(1 < verboseLevel) {
      G4cout << "### SetCrossSectionBiasingFactor: for " 
             << " process " << GetProcessName()
             << " biasFactor= " << f << " weightFlag= " << flag 
             << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ActivateForcedInteraction(G4double length, 
                                                     const G4String& region,
                                                     G4bool flag)
{
  if(nullptr == biasManager) { biasManager = new G4EmBiasingManager(); }
  if(1 < verboseLevel) {
    G4cout << "### ActivateForcedInteraction: for " 
           << " process " << GetProcessName()
           << " length(mm)= " << length/mm
           << " in G4Region <" << region
           << "> weightFlag= " << flag 
           << G4endl; 
  }
  weightFlag = flag;
  biasManager->ActivateForcedInteraction(length, region);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4VEnergyLossProcess::ActivateSecondaryBiasing(const G4String& region, 
                                               G4double factor, 
                                               G4double energyLimit)
{
  if (0.0 <= factor) {
    // Range cut can be applied only for e-
    if(0.0 == factor && secondaryParticle != G4Electron::Electron())
      { return; }

    if(nullptr == biasManager) { biasManager = new G4EmBiasingManager(); }
    biasManager->ActivateSecondaryBiasing(region, factor, energyLimit);
    if(1 < verboseLevel) {
      G4cout << "### ActivateSecondaryBiasing: for " 
             << " process " << GetProcessName()
             << " factor= " << factor
             << " in G4Region <" << region 
             << "> energyLimit(MeV)= " << energyLimit/MeV
             << G4endl; 
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetIonisation(G4bool val)
{
  isIonisation = val;
  aGPILSelection = (val) ? CandidateForSelection : NotCandidateForSelection;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

 void G4VEnergyLossProcess::SetLinearLossLimit(G4double val)
{
  if(0.0 < val && val < 1.0) { 
    linLossLimit = val;
    actLinLossLimit = true; 
  } else { PrintWarning("SetLinearLossLimit", val); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetStepFunction(G4double v1, G4double v2)
{
  if(0.0 < v1 && 0.0 < v2) { 
    dRoverRange = std::min(1.0, v1);
    finalRange = std::min(v2, 1.e+50);
  } else {
    PrintWarning("SetStepFunctionV1", v1); 
    PrintWarning("SetStepFunctionV2", v2); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLowestEnergyLimit(G4double val)
{
  if(1.e-18 < val && val < 1.e+50) { lowestKinEnergy = val; }
  else { PrintWarning("SetLowestEnergyLimit", val); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetDEDXBinning(G4int n)
{
  if(2 < n && n < 1000000000) { 
    nBins = n; 
    actBinning = true;
  } else {
    G4double e = (G4double)n;
    PrintWarning("SetDEDXBinning", e); 
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMinKinEnergy(G4double e)
{
  if(1.e-18 < e && e < maxKinEnergy) { 
    minKinEnergy = e; 
    actMinKinEnergy = true;
  } else { PrintWarning("SetMinKinEnergy", e); } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetMaxKinEnergy(G4double e)
{
  if(minKinEnergy < e && e < 1.e+50) { 
    maxKinEnergy = e;
    actMaxKinEnergy = true;
    if(e < maxKinEnergyCSDA) { maxKinEnergyCSDA = e; }
  } else { PrintWarning("SetMaxKinEnergy", e); } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::PrintWarning(const G4String& tit, G4double val) const
{
  G4String ss = "G4VEnergyLossProcess::" + tit; 
  G4ExceptionDescription ed;
  ed << "Parameter is out of range: " << val 
     << " it will have no effect!\n" << "  Process " 
     << GetProcessName() << "  nbins= " << nBins 
     << " Emin(keV)= " << minKinEnergy/keV 
     << " Emax(GeV)= " << maxKinEnergy/GeV;
  G4Exception(ss, "em0044", JustWarning, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::ProcessDescription(std::ostream& out) const
{
  if(nullptr != particle) { StreamInfo(out, *particle, true); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
