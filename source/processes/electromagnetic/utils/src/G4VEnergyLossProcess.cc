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
#include "G4EmUtility.hh"
#include "G4EmTableUtil.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Electron.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4SafetyHelper.hh"
#include "G4EmDataHandler.hh"
#include "G4TransportationManager.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4VSubCutProducer.hh"
#include "G4EmBiasingManager.hh"
#include "G4Log.hh"
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

namespace 
{
  G4String tnames[7] =
    {"DEDX","Ionisation","DEDXnr","CSDARange","Lambda","Range","InverseRange"};
}


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

  invLambdaFactor = 1.0/lambdaFactor;

  // default linear loss limit
  finalRange = 1.*CLHEP::mm;

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
  isMaster = lManager->IsMaster();

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  theDensityFactor = bld->GetDensityFactors();
  theDensityIdx = bld->GetCoupleIndexes();

  scTracks.reserve(10);
  secParticles.reserve(12);
  emModels = new std::vector<G4VEmModel*>;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEnergyLossProcess::~G4VEnergyLossProcess()
{
  if (isMaster) {
    if(nullptr == baseParticle) { delete theData; }
    delete theEnergyOfCrossSectionMax;
    if(nullptr != fXSpeaks) {
      for(auto const & v : *fXSpeaks) { delete v; }
      delete fXSpeaks;
    }
  }
  delete modelManager;
  delete biasManager;
  delete scoffRegions;
  delete emModels;
  lManager->DeRegister(this);
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
  G4VEmFluctuationModel* afluc = (nullptr == fluc) ? fluctModel : fluc;
  modelManager->AddEmModel(order, ptr, afluc, region);
  ptr->SetParticleChange(pParticleChange, afluc);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetEmModel(G4VEmModel* ptr, G4int)
{
  if(nullptr == ptr) { return; }
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
  particle = G4EmTableUtil::CheckIon(this, &part, particle,
                                     verboseLevel, isIon);

  if( particle != &part ) {
    if(!isIon) { lManager->RegisterExtraParticle(&part, this); }
    if(1 < verboseLevel) {
      G4cout << "### G4VEnergyLossProcess::PreparePhysicsTable()"
             << " interrupted for "
             << part.GetParticleName() << "  isIon=" << isIon << G4endl;
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
  lambdaFactor = theParameters->LambdaFactor();
  invLambdaFactor = 1.0/lambdaFactor;
  if(isMaster) { SetVerboseLevel(theParameters->Verbose()); }
  else { SetVerboseLevel(theParameters->WorkerVerbose()); }
  // integral option may be disabled
  if(!theParameters->Integral()) { fXSType = fEmNoIntegral; }

  theParameters->DefineRegParamForLoss(this);

  fRangeEnergy = 0.0;

  G4double initialCharge = particle->GetPDGCharge();
  G4double initialMass   = particle->GetPDGMass();

  theParameters->FillStepFunction(particle, this);

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
    if(nullptr == theData) { theData = new G4EmDataHandler(7); }

    if(nullptr != theDEDXTable && isIonisation) {
      if(nullptr != theIonisationTable && theDEDXTable != theIonisationTable) {
	theData->CleanTable(0);
	theDEDXTable = theIonisationTable;
	theIonisationTable = nullptr;
      }
    }
    
    theDEDXTable = theData->MakeTable(theDEDXTable, 0);
    bld->InitialiseBaseMaterials(theDEDXTable);
    theData->UpdateTable(theIonisationTable, 1);

    if (theParameters->BuildCSDARange()) {
      theDEDXunRestrictedTable = theData->MakeTable(2);
      if(isIonisation) { theCSDARangeTable = theData->MakeTable(3); }
    }

    theLambdaTable = theData->MakeTable(4);
    if(isIonisation) {
      theRangeTableForLoss = theData->MakeTable(5);
      theInverseRangeTable = theData->MakeTable(6);
    }
  }

  // forced biasing
  if(nullptr != biasManager) { 
    biasManager->Initialise(part,GetProcessName(),verboseLevel); 
    biasFlag = false; 
  }
  baseMat = bld->GetBaseMaterialFlag();
  numberOfModels = modelManager->NumberOfModels();
  currentModel = modelManager->GetModel(0);
  G4EmTableUtil::UpdateModels(this, modelManager, maxKinEnergy,
                              numberOfModels, secID, biasID,
                              mainSecondaries, baseMat, isMaster,
                              theParameters->UseAngularGeneratorForIonisation());
  theCuts = modelManager->Initialise(particle, secondaryParticle,
                                     verboseLevel);
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
      const auto masterProcess =
        static_cast<const G4VEnergyLossProcess*>(GetMasterProcess());

      numberOfModels = modelManager->NumberOfModels();
      G4EmTableUtil::BuildLocalElossProcess(this, masterProcess,
                                            particle, numberOfModels);
      tablesAreBuilt = true;  
      baseMat = masterProcess->UseBaseMaterial();
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
                           num == "GenericIon"|| num == "alpha+" ))) { 
    StreamInfo(G4cout, part); 
  }
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
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossProcess::BuildDEDXTable() of type " << tType
           << " for " << GetProcessName()
           << " and " << particle->GetParticleName() << G4endl;
  }
  if(nullptr == table) { return table; }

  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  G4EmTableUtil::BuildDEDXTable(this, particle, modelManager, bld,
                                table, minKinEnergy, emax, bin,
                                verboseLevel, tType, spline);
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEnergyLossProcess::BuildLambdaTable(G4EmTableType)
{
  if(nullptr == theLambdaTable) { return theLambdaTable; }

  G4double scale = theParameters->MaxKinEnergy()/theParameters->MinKinEnergy();
  G4int nbin = 
    theParameters->NumberOfBinsPerDecade()*G4lrint(std::log10(scale));
  scale = nbin/G4Log(scale);
  
  G4LossTableBuilder* bld = lManager->GetTableBuilder();
  G4EmTableUtil::BuildLambdaTable(this, particle, modelManager,
                                  bld, theLambdaTable, theCuts,
                                  minKinEnergy, maxKinEnergy, scale,
                                  verboseLevel, spline);
  return theLambdaTable;
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
    for(std::size_t i=0; i<7; ++i) {
      auto ta = theData->Table(i);
      out << "      " << tnames[i] << " address: " << ta << G4endl; 
      if(nullptr != ta) { out << *ta << G4endl; }
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
  // reset parameters for the new track
  theNumberOfInteractionLengthLeft = -1.0;
  mfpKinEnergy = DBL_MAX;
  currentCouple = nullptr;

  // reset ion
  if(isIon) {
    const G4double newmass = track->GetDefinition()->GetPDGMass();
    if(nullptr != baseParticle) {
      massRatio = baseParticle->GetPDGMass()/newmass;
      logMassRatio = G4Log(massRatio);
    } else if(isIon) {
      massRatio = CLHEP::proton_mass_c2/newmass;
      logMassRatio = G4Log(massRatio);
    } else {
      massRatio = 1.0;
      logMassRatio = 0.0;
    }
  }  
  // forced biasing only for primary particles
  if(nullptr != biasManager) {
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
    }
    if (lossFluctuationFlag) {
      auto fluc = currentModel->GetModelOfFluctuations();
      fluc->SetParticleAndCharge(track.GetDefinition(), q2);
    }
  }

  // forced biasing only for primary particles
  if(biasManager) {
    if(0 == track.GetParentID() && biasFlag && 
       biasManager->ForcedInteractionRegion((G4int)currentCoupleIndex)) {
      return biasManager->GetStepLimit((G4int)currentCoupleIndex, previousStepSize);
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
  if (verboseLevel>2) {
    G4cout << "G4VEnergyLossProcess::PostStepGetPhysicalInteractionLength ";
    G4cout << "[ " << GetProcessName() << "]" << G4endl; 
    G4cout << " for " << track.GetDefinition()->GetParticleName() 
           << " in Material  " <<  currentMaterial->GetName()
           << " Ekin(MeV)= " << preStepKinEnergy/MeV 
           << " track material: " << track.GetMaterial()->GetName()
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
    if(e*invLambdaFactor < mfpKinEnergy) {
      mfpKinEnergy = e;
      preStepLambda = GetLambdaForScaledEnergy(e, loge); 
    }

    // cross section has one peak
  } else if(fXSType == fEmOnePeak) {
    const G4double epeak = (*theEnergyOfCrossSectionMax)[basedCoupleIndex];
    if(e <= epeak) {
      if(e*invLambdaFactor < mfpKinEnergy) {
        mfpKinEnergy = e;
        preStepLambda = GetLambdaForScaledEnergy(e, loge); 
      }
    } else if(e < mfpKinEnergy) { 
      const G4double e1 = std::max(epeak, e*lambdaFactor);
      mfpKinEnergy = e1;
      preStepLambda = GetLambdaForScaledEnergy(e1);
    }

    // cross section has more than one peaks
  } else if(fXSType == fEmTwoPeaks) {
    G4TwoPeaksXS* xs = (*fXSpeaks)[basedCoupleIndex];
    const G4double e1peak = xs->e1peak;

    // below the 1st peak
    if(e <= e1peak) {
      if(e*invLambdaFactor < mfpKinEnergy) {
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
        mfpKinEnergy = e1;
        preStepLambda = GetLambdaForScaledEnergy(e1); 
      }
      return;
    }
    const G4double e2peak = xs->e2peak;
    // above the deep, below 2nd peak
    if(e <= e2peak) {
      if(e*invLambdaFactor < mfpKinEnergy) {
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
        mfpKinEnergy = e1;
        preStepLambda = GetLambdaForScaledEnergy(e1); 
      }
      return;
    }
    const G4double e3peak = xs->e3peak;
    // above the deep, below 3d peak
    if(e <= e3peak) {
      if(e*invLambdaFactor < mfpKinEnergy) {
        mfpKinEnergy = e;
        preStepLambda = GetLambdaForScaledEnergy(e, loge); 
      }
      return;
    }
    // above 3d peak
    if(e <= mfpKinEnergy) { 
      const G4double e1 = std::max(e3peak, e*lambdaFactor);
      mfpKinEnergy = e1;
      preStepLambda = GetLambdaForScaledEnergy(e1); 
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
           << GetProcessName() << " and particle " << d->GetParticleName()
           << "  eScaled(MeV)=" << preStepScaledEnergy/MeV
           << "  range(mm)=" << fRange/mm << "  s(mm)=" << length/mm
           << "  rf=" << reduceFactor << "  q^2=" << chargeSqRatio
           << " md=" << d->GetPDGMass() << "  status=" << track.GetTrackStatus()
           << "  " << track.GetMaterial()->GetName() << G4endl;
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
                                              eloss, (G4int)currentCoupleIndex);
      if(scTracks.size() > 0) { FillSecondariesAlongStep(weight); }
      eloss = std::max(eloss, 0.0);
    }
    fParticleChange.SetProposedKineticEnergy(0.0);
    fParticleChange.ProposeLocalEnergyDeposit(eloss);
    return &fParticleChange;
  }
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
             << " eloss(MeV)= " << eloss/MeV << " eloss0(MeV)= "
             << GetDEDXForScaledEnergy(preStepScaledEnergy)*length/MeV
             << " lim(MeV)= " << preStepKinEnergy*linLossLimit/MeV
             << G4endl;
    */
  }

  /*
  if(-1 < verboseLevel ) {
    G4cout << "Before fluct: eloss(MeV)= " << eloss/MeV
           << " e-eloss= " << preStepKinEnergy-eloss
           << " step(mm)= " << length/mm << " range(mm)= " << fRange/mm
           << " fluct= " << lossFluctuationFlag << G4endl;
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
             << " massRatio= " << massRatio << " tmax= " << tmax << G4endl;
    */
  }

  // deexcitation
  if (useDeexcitation) {
    G4double esecfluo = preStepKinEnergy;
    G4double de = esecfluo;
    atomDeexcitation->AlongStepDeexcitation(scTracks, step, 
                                            de, (G4int)currentCoupleIndex);

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
  const std::size_t n0 = scTracks.size();
  G4double weight = wt;
  // weight may be changed by biasing manager
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion((G4int)currentCoupleIndex)) {
      weight *=
        biasManager->ApplySecondaryBiasing(scTracks, (G4int)currentCoupleIndex);
    }
  } 

  // fill secondaries
  const std::size_t n = scTracks.size();
  fParticleChange.SetNumberOfSecondaries((G4int)n);

  for(std::size_t i=0; i<n; ++i) {
    G4Track* t = scTracks[i];
    if(nullptr != t) {
      t->SetWeight(weight); 
      pParticleChange->AddSecondary(t);
      if(i >= n0) { t->SetCreatorModelID(biasID); }
    }
  }
  scTracks.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEnergyLossProcess::PostStepDoIt(const G4Track& track,
                                                      const G4Step& step)
{
  // clear number of interaction lengths in any case
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
  if(1 < verboseLevel) {
    G4cout<<GetProcessName()<<" PostStepDoIt: E(MeV)= "<< finalT/MeV<< G4endl;
  }
  */
  // forced process - should happen only once per track
  if(biasFlag) {
    if(biasManager->ForcedInteractionRegion((G4int)currentCoupleIndex)) {
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

    // if both lg and lx are zero then no interaction
    if(preStepLambda*G4UniformRand() >= lx) {
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
  currentModel->SampleSecondaries(&secParticles, currentCouple, dp, tcut);

  const G4int num0 = (G4int)secParticles.size();

  // bremsstrahlung splitting or Russian roulette  
  if(biasManager) {
    if(biasManager->SecondaryBiasingRegion((G4int)currentCoupleIndex)) {
      G4double eloss = 0.0;
      weight *= biasManager->ApplySecondaryBiasing(
                                      secParticles,
                                      track, currentModel, 
                                      &fParticleChange, eloss,
                                      (G4int)currentCoupleIndex, tcut, 
                                      step.GetPostStepPoint()->GetSafety());
      if(eloss > 0.0) {
        eloss += fParticleChange.GetLocalEnergyDeposit();
        fParticleChange.ProposeLocalEnergyDeposit(eloss);
      }
    }
  }

  // save secondaries
  const G4int num = (G4int)secParticles.size();
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
       const G4ParticleDefinition* part, const G4String& dir, G4bool ascii)
{
  if (!isMaster || nullptr != baseParticle || part != particle ) return true;
  for(std::size_t i=0; i<7; ++i) {
    if(nullptr != theData->Table(i)) {
      if(1 < verboseLevel) {
	G4cout << "G4VEnergyLossProcess::StorePhysicsTable i=" << i
	       << "  " << particle->GetParticleName()
	       << "  " << GetProcessName()
	       << "  " << tnames[i] << "  " << theData->Table(i) << G4endl;
      }
      if(!G4EmTableUtil::StoreTable(this, part, theData->Table(i), 
				    dir, tnames[i], verboseLevel, ascii)) { 
	return false;
      }
    }
  }
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4bool 
G4VEnergyLossProcess::RetrievePhysicsTable(const G4ParticleDefinition* part, 
                                           const G4String& dir, G4bool ascii)
{
  if (!isMaster || nullptr != baseParticle || part != particle ) return true;
  for(std::size_t i=0; i<7; ++i) {
    if(!G4EmTableUtil::RetrieveTable(this, part, theData->Table(i), dir, tnames[i],
                                     verboseLevel, ascii, spline)) { 
      return false;
    }
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
  if(1 < verboseLevel) {
    G4cout << "### Set DEDX table " << p << "  " << theDEDXTable
	   << "  " <<  theDEDXunRestrictedTable << "  " << theIonisationTable
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() 
	   << " type=" << tType << " isIonisation:" << isIonisation << G4endl;
  }
  if(fTotal == tType) {
    theDEDXunRestrictedTable = p;
  } else if(fRestricted == tType) {
    theDEDXTable = p;
    if(isMaster && nullptr == baseParticle) {
      theData->UpdateTable(theDEDXTable, 0);
    }
  } else if(fIsIonisation == tType) {
    theIonisationTable = p;
    if(isMaster && nullptr == baseParticle) {
      theData->UpdateTable(theIonisationTable, 1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetCSDARangeTable(G4PhysicsTable* p)
{
  theCSDARangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetRangeTableForLoss(G4PhysicsTable* p)
{
  theRangeTableForLoss = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetInverseRangeTable(G4PhysicsTable* p)
{
  theInverseRangeTable = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetLambdaTable(G4PhysicsTable* p)
{
  if(1 < verboseLevel) {
    G4cout << "### Set Lambda table " << p << " " << theLambdaTable 
           << " for " << particle->GetParticleName()
           << " and process " << GetProcessName() << G4endl;
  }
  theLambdaTable = p;
  tablesAreBuilt = true;

  if(isMaster && nullptr != p) {
    delete theEnergyOfCrossSectionMax;
    theEnergyOfCrossSectionMax = nullptr;
    if(fEmTwoPeaks == fXSType) {
      if(nullptr != fXSpeaks) { 
	for(auto & ptr : *fXSpeaks) { delete ptr; }
	delete fXSpeaks;
      }
      G4LossTableBuilder* bld = lManager->GetTableBuilder();
      fXSpeaks = G4EmUtility::FillPeaksStructure(p, bld);
      if(nullptr == fXSpeaks) { fXSType = fEmOnePeak; }
    }
    if(fXSType == fEmOnePeak) { 
      theEnergyOfCrossSectionMax = G4EmUtility::FindCrossSectionMax(p);
      if(nullptr == theEnergyOfCrossSectionMax) { fXSType = fEmIncreasing; }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetEnergyOfCrossSectionMax(std::vector<G4double>* p)
{
  theEnergyOfCrossSectionMax = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEnergyLossProcess::SetTwoPeaksXS(std::vector<G4TwoPeaksXS*>* ptr)
{
  fXSpeaks = ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4Element* G4VEnergyLossProcess::GetCurrentElement() const
{
  return (nullptr != currentModel) 
    ? currentModel->GetCurrentElement(currentMaterial) : nullptr;
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
