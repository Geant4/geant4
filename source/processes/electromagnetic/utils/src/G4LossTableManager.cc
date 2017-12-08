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
// $Id: G4LossTableManager.cc 107364 2017-11-09 10:53:25Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4LossTableManager
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 03.01.2002
//
// Modifications:
//
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 15-02-03 Lambda table can be scaled (V.Ivanchenko)
// 17-02-03 Fix problem of store/restore tables (V.Ivanchenko)
// 10-03-03 Add Ion registration (V.Ivanchenko)
// 25-03-03 Add deregistration (V.Ivanchenko)
// 02-04-03 Change messenger (V.Ivanchenko)
// 26-04-03 Fix retrieve tables (V.Ivanchenko)
// 13-05-03 Add calculation of precise range (V.Ivanchenko)
// 23-07-03 Add exchange with G4EnergyLossTables (V.Ivanchenko)
// 05-10-03 Add G4VEmProcesses registration and Verbose command (V.Ivanchenko)
// 17-10-03 Add SetParameters method (V.Ivanchenko)
// 23-10-03 Add control on inactive processes (V.Ivanchenko)
// 04-11-03 Add checks in RetrievePhysicsTable (V.Ivanchenko)
// 12-11-03 G4EnergyLossSTD -> G4EnergyLossProcess (V.Ivanchenko)
// 14-01-04 Activate precise range calculation (V.Ivanchenko)
// 10-03-04 Fix a problem of Precise Range table (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 13-01-04 Fix problem which takes place for inactivate eIoni (V.Ivanchenko)
// 25-01-04 Fix initialisation problem for ions (V.Ivanchenko)
// 11-03-05 Shift verbose level by 1 (V.Ivantchenko)
// 10-01-06 PreciseRange -> CSDARange (V.Ivantchenko)
// 20-01-06 Introduce G4EmTableType to remove repeating code (VI)
// 23-03-06 Set flag isIonisation (VI)
// 10-05-06 Add methods  SetMscStepLimitation, FacRange and MscFlag (VI)
// 22-05-06 Add methods  Set/Get bremsTh (VI)
// 05-06-06 Do not clear loss_table map between runs (VI)
// 16-01-07 Create new energy loss table for e+,e-,mu+,mu- and 
//          left ionisation table for further usage (VI)
// 12-02-07 Add SetSkin, SetLinearLossLimit (V.Ivanchenko)
// 18-06-07 Move definition of msc parameters to G4EmProcessOptions (V.Ivanchenko)
// 21-02-08 Added G4EmSaturation (V.Ivanchenko)
// 12-04-10 Added PreparePhysicsTables and BuildPhysicsTables entries (V.Ivanchenko)
// 04-06-13 (V.Ivanchenko) Adaptation for MT mode; new method LocalPhysicsTables; 
//          ions expect G4GenericIon are not included in the map of energy loss
//          processes for performnc reasons  
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LossTableManager.hh"
#include "G4SystemOfUnits.hh"

#include "G4VMultipleScattering.hh"
#include "G4VEmProcess.hh"

#include "G4EmParameters.hh"
#include "G4EmSaturation.hh"
#include "G4EmConfigurator.hh"
#include "G4ElectronIonPair.hh"

#include "G4PhysicsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4EmTableType.hh"
#include "G4Region.hh"
#include "G4PhysicalConstants.hh"
#include "G4Threading.hh"

#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4OpticalPhoton.hh"
#include "G4Neutron.hh"
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4ThreadLocal G4LossTableManager* G4LossTableManager::instance = nullptr;

G4LossTableManager* G4LossTableManager::Instance()
{
  if(!instance) {
    static G4ThreadLocalSingleton<G4LossTableManager> inst;
    instance = inst.Instance();
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::~G4LossTableManager()
{
  //G4cout << "### G4LossTableManager::~G4LossTableManager() "<< this << G4endl;
  for (G4int i=0; i<n_loss; ++i) {
    //G4cout << "### eloss #" << i << G4endl;
    if( loss_vector[i] ) { 
      delete loss_vector[i]; 
    }
  }
  size_t msc = msc_vector.size();
  for (size_t j=0; j<msc; ++j) {
    if( msc_vector[j] ) { delete msc_vector[j]; }
  }
  size_t emp = emp_vector.size();
  for (size_t k=0; k<emp; ++k) {
    if( emp_vector[k] ) { delete emp_vector[k]; }
  }
  size_t mod = mod_vector.size();
  size_t fmod = fmod_vector.size();
  //G4cout << "   Nmod" << mod << "  Nfluc= " << fmod << G4endl;
  for (size_t a=0; a<mod; ++a) {
    //G4cout << "Delete model #" << a << "  " <<  mod_vector[a] << G4endl;
    if( nullptr != mod_vector[a] ) { 
      for (size_t b=0; b<fmod; ++b) {
        if((G4VEmModel*)(fmod_vector[b]) == mod_vector[a]) {
          fmod_vector[b] = nullptr;
        }
      }
      delete mod_vector[a]; 
      mod_vector[a] = nullptr;
    }
  }
  for (size_t b=0; b<fmod; ++b) {
    if( fmod_vector[b] ) { delete fmod_vector[b]; }
  }
  Clear();
  delete tableBuilder;
  delete emCorrections;
  delete emConfigurator;
  delete emElectronIonPair;
  delete atomDeexcitation;
  delete subcutProducer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::G4LossTableManager()
{
  //G4cout << "### G4LossTableManager::G4LossTableManager() " << this << G4endl;
  theParameters = G4EmParameters::Instance();
  n_loss = 0;
  run = -1;
  startInitialisation = false;
  all_tables_are_built = false;
  currentLoss = nullptr;
  currentParticle = nullptr;
  firstParticle = nullptr;
  isMaster = true;
  verbose = theParameters->Verbose();
  theElectron  = G4Electron::Electron();
  theGenericIon= nullptr;
  if(G4Threading::IsWorkerThread()) { 
    verbose = theParameters->WorkerVerbose();
    isMaster = false;
  }  
  tableBuilder = new G4LossTableBuilder();
  emCorrections= new G4EmCorrections(verbose);
  emConfigurator = nullptr;
  emElectronIonPair = nullptr;
  atomDeexcitation = nullptr;
  subcutProducer = nullptr;
  gammaShark = nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Clear()
{
  all_tables_are_built = false;
  currentLoss = nullptr;
  currentParticle = nullptr;
  if(n_loss)
    {
      dedx_vector.clear();
      range_vector.clear();
      inv_range_vector.clear();
      loss_map.clear();
      loss_vector.clear();
      part_vector.clear();
      base_part_vector.clear();
      tables_are_built.clear();
      isActive.clear();
      n_loss = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEnergyLossProcess* p)
{
  if(!p) { return; }
  for (G4int i=0; i<n_loss; ++i) {
    if(loss_vector[i] == p) { return; }
  }
  if(verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEnergyLossProcess : " 
           << p->GetProcessName() << "  idx= " << n_loss << G4endl;
  }
  ++n_loss;
  loss_vector.push_back(p);
  part_vector.push_back(nullptr);
  base_part_vector.push_back(nullptr);
  dedx_vector.push_back(nullptr);
  range_vector.push_back(nullptr);
  inv_range_vector.push_back(nullptr);
  tables_are_built.push_back(false);
  isActive.push_back(true);
  all_tables_are_built = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::ResetParameters()
{
  verbose = theParameters->Verbose();
  if(!isMaster) {
    verbose = theParameters->WorkerVerbose();
  }
  tableBuilder->SetSplineFlag(theParameters->Spline());
  tableBuilder->SetInitialisationFlag(false); 
  emCorrections->SetVerbose(verbose); 
  if(emConfigurator) { emConfigurator->SetVerbose(verbose); };
  if(emElectronIonPair) { emElectronIonPair->SetVerbose(verbose); };
  if(atomDeexcitation) {
    atomDeexcitation->SetVerboseLevel(verbose);
    atomDeexcitation->InitialiseAtomicDeexcitation();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEnergyLossProcess* p)
{
  if(!p) { return; }
  for (G4int i=0; i<n_loss; ++i) {
    if(loss_vector[i] == p) { loss_vector[i] = nullptr; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VMultipleScattering* p)
{
  if(!p) { return; }
  G4int n = msc_vector.size();
  for (G4int i=0; i<n; ++i) {
    if(msc_vector[i] == p) { return; }
  }
  if(verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VMultipleScattering : " 
           << p->GetProcessName() << "  idx= " << msc_vector.size() << G4endl;
  }
  msc_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VMultipleScattering* p)
{
  if(!p) { return; }
  size_t msc = msc_vector.size();
  for (size_t i=0; i<msc; ++i) {
    if(msc_vector[i] == p) { msc_vector[i] = nullptr; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmProcess* p)
{
  if(!p) { return; }
  G4int n = emp_vector.size();
  for (G4int i=0; i<n; ++i) {
    if(emp_vector[i] == p) { return; }
  }
  if(verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEmProcess : " 
           << p->GetProcessName() << "  idx= " << emp_vector.size() << G4endl;
  }
  emp_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmProcess* p)
{
  if(!p) { return; }
  size_t emp = emp_vector.size();
  for (size_t i=0; i<emp; ++i) {
    if(emp_vector[i] == p) { emp_vector[i] = nullptr; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmModel* p)
{
  mod_vector.push_back(p);
  if(verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEmModel : " 
           << p->GetName() << "  " << p << "  " << mod_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmModel* p)
{
  //G4cout << "G4LossTableManager::DeRegister G4VEmModel : " << p << G4endl;
  size_t n = mod_vector.size();
  for (size_t i=0; i<n; ++i) {
    if(mod_vector[i] == p) { 
      mod_vector[i] = nullptr; 
      break;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmFluctuationModel* p)
{
  fmod_vector.push_back(p);
  if(verbose > 1) {
    G4cout << "G4LossTableManager::Register G4VEmFluctuationModel : " 
           << p->GetName() << "  " << fmod_vector.size() << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmFluctuationModel* p)
{
  size_t n = fmod_vector.size();
  for (size_t i=0; i<n; ++i) {
    if(fmod_vector[i] == p) { fmod_vector[i] = nullptr; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::RegisterExtraParticle(
     const G4ParticleDefinition* part,
     G4VEnergyLossProcess* p)
{ 
  if(!p || !part) { return; }
  for (G4int i=0; i<n_loss; ++i) {
    if(loss_vector[i] == p) { return; }
  }
  if(verbose > 1) {
    G4cout << "G4LossTableManager::RegisterExtraParticle "
           << part->GetParticleName() << "  G4VEnergyLossProcess : " 
           << p->GetProcessName() << "  idx= " << n_loss << G4endl;
  }
  ++n_loss;
  loss_vector.push_back(p);
  part_vector.push_back(part);
  base_part_vector.push_back(p->BaseParticle());
  dedx_vector.push_back(nullptr);
  range_vector.push_back(nullptr);
  inv_range_vector.push_back(nullptr);
  tables_are_built.push_back(false);
  all_tables_are_built = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* particle,
                                        G4VEnergyLossProcess* p,
                                        G4bool theMaster)
{
  if (1 < verbose) {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " 
           << particle->GetParticleName() 
           << " and " << p->GetProcessName() << " run= " << run 
           << "   loss_vector " << loss_vector.size() << G4endl;
  }

  isMaster = theMaster;

  if(!startInitialisation) { 
    ResetParameters();
    if (1 < verbose) {
      G4cout << "====== G4LossTableManager::PreparePhysicsTable start ====="
             << G4endl;
    } 
  }

  // start initialisation for the first run
  if( -1 == run ) {
    if(emConfigurator) { emConfigurator->PrepareModels(particle, p); }

    // initialise particles for given process
    for (G4int j=0; j<n_loss; ++j) {
      if (p == loss_vector[j] && !part_vector[j]) { 
        part_vector[j] = particle;
        if(particle->GetParticleName() == "GenericIon") {
          theGenericIon = particle;
        }
      }
    }
  }
  startInitialisation = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* particle,
                                        G4VEmProcess* p, G4bool theMaster)
{
  if (1 < verbose) {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " 
           << particle->GetParticleName() 
           << " and " << p->GetProcessName() << G4endl;
  }
  isMaster = theMaster;

  if(!startInitialisation) { 
    ResetParameters(); 
    if (1 < verbose) {
      G4cout << "====== G4LossTableManager::PreparePhysicsTable start ====="
             << G4endl;
    } 
  }

  // start initialisation for the first run
  if( -1 == run ) {
    if(emConfigurator) { emConfigurator->PrepareModels(particle, p); }
  }
  startInitialisation = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4LossTableManager::PreparePhysicsTable(const G4ParticleDefinition* particle,
                                        G4VMultipleScattering* p,
                                        G4bool theMaster)
{
  if (1 < verbose) {
    G4cout << "G4LossTableManager::PreparePhysicsTable for " 
           << particle->GetParticleName() 
           << " and " << p->GetProcessName() << G4endl;
  }

  isMaster = theMaster;

  if(!startInitialisation) { 
    ResetParameters();
    if (1 < verbose) {
      G4cout << "====== G4LossTableManager::PreparePhysicsTable start ====="
             << G4endl;
    } 
  }

  // start initialisation for the first run
  if( -1 == run ) {
    if(emConfigurator) { emConfigurator->PrepareModels(particle, p); }
  } 
  startInitialisation = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void 
G4LossTableManager::BuildPhysicsTable(const G4ParticleDefinition*)
{
  if(-1 == run && startInitialisation) {
    if(emConfigurator) { emConfigurator->Clear(); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::LocalPhysicsTables(
     const G4ParticleDefinition* aParticle,
     G4VEnergyLossProcess* p)
{
  if(1 < verbose) {
    G4cout << "### G4LossTableManager::LocalPhysicsTable() for "
           << aParticle->GetParticleName()
	   << " and process " << p->GetProcessName()
           << G4endl;
  }

  if(-1 == run && startInitialisation) {
    if(emConfigurator) { emConfigurator->Clear(); }
    firstParticle = aParticle; 
  }

  if(startInitialisation) {
    ++run;
    if(1 < verbose) {
      G4cout << "===== G4LossTableManager::LocalPhysicsTable() for run "
             << run << " =====" << G4endl;
    }
    currentParticle = nullptr;
    startInitialisation = false;
    for (G4int i=0; i<n_loss; ++i) {
      if(loss_vector[i]) {
        tables_are_built[i] = false;
      } else {
        tables_are_built[i] = true;
        part_vector[i] = nullptr;
      }
    }
  }

  all_tables_are_built= true;
  for (G4int i=0; i<n_loss; ++i) {
    if(p == loss_vector[i]) {
      tables_are_built[i] = true;
      isActive[i] = true;
      part_vector[i] = p->Particle(); 
      base_part_vector[i] = p->BaseParticle(); 
      dedx_vector[i] = p->DEDXTable();
      range_vector[i] = p->RangeTableForLoss();
      inv_range_vector[i] = p->InverseRangeTable();
      if(0 == run && p->IsIonisationProcess()) {
        loss_map[part_vector[i]] = p;
	//G4cout << "G4LossTableManager::LocalPhysicsTable " << part_vector[i]->GetParticleName()
	//     << " added to map " << p << G4endl;
      }

      if(1 < verbose) { 
        G4cout << i <<".   "<< p->GetProcessName(); 
        if(part_vector[i]) {
          G4cout << "  for "  << part_vector[i]->GetParticleName();
        }
        G4cout << "  active= " << isActive[i]
               << "  table= " << tables_are_built[i]
               << "  isIonisation= " << p->IsIonisationProcess()
               << G4endl;
      }
      break;
    } else if(!tables_are_built[i]) {
      all_tables_are_built = false;
    }
  }

  if(1 < verbose) {
    G4cout << "### G4LossTableManager::LocalPhysicsTable end"
           << G4endl;
  }
  if(all_tables_are_built) { 
    if(1 < verbose) {
      G4cout << "%%%%% All dEdx and Range tables for worker are ready for run "
             << run << " %%%%%" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::BuildPhysicsTable(
     const G4ParticleDefinition* aParticle,
     G4VEnergyLossProcess* p)
{
  if(1 < verbose) {
    G4cout << "### G4LossTableManager::BuildPhysicsTable() for "
           << aParticle->GetParticleName()
           << " and process " << p->GetProcessName() << G4endl;
  }
  // clear configurator
  if(-1 == run && startInitialisation) {
    if(emConfigurator) { emConfigurator->Clear(); }
    firstParticle = aParticle; 
  }
  if(startInitialisation) {
    ++run;
    if(1 < verbose) {
      G4cout << "===== G4LossTableManager::BuildPhysicsTable() for run "
             << run << " ===== " << atomDeexcitation << G4endl;
    }
    currentParticle = nullptr;
    all_tables_are_built= true;
  }

  // initialisation before any table is built
  if ( startInitialisation && aParticle == firstParticle ) {

    startInitialisation = false;
    if(1 < verbose) {
      G4cout << "### G4LossTableManager start initilisation for first particle "
             << firstParticle->GetParticleName() 
             << G4endl;
    }
    for (G4int i=0; i<n_loss; ++i) {
      G4VEnergyLossProcess* el = loss_vector[i];

      if(el) {
        isActive[i] = true;
        base_part_vector[i] = el->BaseParticle(); 
        tables_are_built[i] = false;
        all_tables_are_built= false;
        if(!isActive[i]) { 
          el->SetIonisation(false); 
          tables_are_built[i] = true;
        }
  
        if(1 < verbose) { 
          G4cout << i <<".   "<< el->GetProcessName();
          if(el->Particle()) {
            G4cout << "  for "  << el->Particle()->GetParticleName();
          }
          G4cout << "  active= " << isActive[i]
                 << "  table= " << tables_are_built[i]
                 << "  isIonisation= " << el->IsIonisationProcess();
          if(base_part_vector[i]) { 
            G4cout << "  base particle " 
                   << base_part_vector[i]->GetParticleName();
          }
          G4cout << G4endl;
        }
      } else {
        tables_are_built[i] = true;
        part_vector[i] = nullptr;
        isActive[i] = false;
      }
    }
  }

 if (all_tables_are_built) { return; }

  // Build tables for given particle
  all_tables_are_built = true;

  for(G4int i=0; i<n_loss; ++i) {
    if(p == loss_vector[i] && !tables_are_built[i] && !base_part_vector[i]) {
      const G4ParticleDefinition* curr_part = part_vector[i];
      if(1 < verbose) {
        G4cout << "### Build Table for " << p->GetProcessName()
               << " and " << curr_part->GetParticleName() 
	       << "  " << tables_are_built[i] << "  " << base_part_vector[i] << G4endl;
      }
      G4VEnergyLossProcess* curr_proc = BuildTables(curr_part);
      if(curr_proc) { 
	CopyTables(curr_part, curr_proc); 
	if(p == curr_proc && 0 == run && p->IsIonisationProcess()) { 
	  loss_map[aParticle] = p; 
	  //G4cout << "G4LossTableManager::BuildPhysicsTable: " << aParticle->GetParticleName()
	  //	 << " added to map " << p <<  G4endl;
	}
      }
    }
    if ( !tables_are_built[i] ) { all_tables_are_built = false; }
  }
  if(1 < verbose) {
    G4cout << "### G4LossTableManager::BuildPhysicsTable end: "
           << "all_tables_are_built= " << all_tables_are_built << " "
           << aParticle->GetParticleName() << " proc: " << p << G4endl;
  }
  if(all_tables_are_built) { 
    if(1 < verbose) {
      G4cout << "%%%%% All dEdx and Range tables are built for master run= " 
             << run << " %%%%%" << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::CopyTables(const G4ParticleDefinition* part,
                                    G4VEnergyLossProcess* base_proc)
{
  for (G4int j=0; j<n_loss; ++j) {

    G4VEnergyLossProcess* proc = loss_vector[j];

    if (!tables_are_built[j] && part == base_part_vector[j]) {
      tables_are_built[j] = true;
      proc->SetDEDXTable(base_proc->IonisationTable(),fRestricted);
      proc->SetDEDXTable(base_proc->DEDXTableForSubsec(),fSubRestricted);
      proc->SetDEDXTable(base_proc->DEDXunRestrictedTable(),fTotal);
      proc->SetCSDARangeTable(base_proc->CSDARangeTable());
      proc->SetRangeTableForLoss(base_proc->RangeTableForLoss());
      proc->SetInverseRangeTable(base_proc->InverseRangeTable());
      proc->SetLambdaTable(base_proc->LambdaTable());
      proc->SetSubLambdaTable(base_proc->SubLambdaTable());
      proc->SetIonisation(base_proc->IsIonisationProcess());
      if(proc->IsIonisationProcess()) { 
        range_vector[j] = base_proc->RangeTableForLoss();
        inv_range_vector[j] = base_proc->InverseRangeTable();
        loss_map[part_vector[j]] = proc; 
	//G4cout << "G4LossTableManager::CopyTable " << part_vector[j]->GetParticleName()
	//     << " added to map " << proc << G4endl;
      }
      if (1 < verbose) {
         G4cout << "For " << proc->GetProcessName()
                << " for " << part_vector[j]->GetParticleName()
                << " base_part= " << part->GetParticleName()
                << " tables are assigned"
                << G4endl;
      }
    }

    if (theElectron == part && theElectron == proc->SecondaryParticle() ) {
      proc->SetSecondaryRangeTable(base_proc->RangeTableForLoss());
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4VEnergyLossProcess* G4LossTableManager::BuildTables(
                      const G4ParticleDefinition* aParticle)
{
  if(1 < verbose) {
    G4cout << "G4LossTableManager::BuildTables() for "
           << aParticle->GetParticleName() << G4endl;
  }

  std::vector<G4PhysicsTable*> t_list;  
  std::vector<G4VEnergyLossProcess*> loss_list;
  std::vector<G4bool> build_flags;
  G4VEnergyLossProcess* em = nullptr;
  G4VEnergyLossProcess* p = nullptr;
  G4int iem = 0;
  G4PhysicsTable* dedx = 0;
  G4int i;

  G4ProcessVector* pvec = 
    aParticle->GetProcessManager()->GetProcessList();
  G4int nvec = pvec->size();

  for (i=0; i<n_loss; ++i) {
    p = loss_vector[i];
    if (p) {
      G4bool yes = (aParticle == part_vector[i]);

      // possible case of process sharing between particle/anti-particle
      if(!yes) {
        G4VProcess* ptr = static_cast<G4VProcess*>(p);
        for(G4int j=0; j<nvec; ++j) {
          //G4cout << "j= " << j << " " << (*pvec)[j] << " " << ptr << G4endl;
          if(ptr == (*pvec)[j]) {
            yes = true;
            break;
          }
        }
      }      
      // process belong to this particle
      if(yes && isActive[i]) {
        if (p->IsIonisationProcess() || !em) {
          em = p;
          iem= i;
        }
        // tables may be shared between particle/anti-particle
        G4bool val = false;
        if (!tables_are_built[i]) {
          val = true;
          dedx = p->BuildDEDXTable(fRestricted);
          //G4cout << "Build DEDX table for " << p->GetProcessName()
          // << " idx= " << i << dedx << " " << dedx->length() << G4endl;
          p->SetDEDXTable(dedx,fRestricted);
          tables_are_built[i] = true;
        } else {
          dedx = p->DEDXTable();
        }
        t_list.push_back(dedx);
        loss_list.push_back(p);
        build_flags.push_back(val);
      }
    }
  }

  G4int n_dedx = t_list.size();
  if (0 == n_dedx || !em) {
    G4cout << "G4LossTableManager WARNING: no DEDX processes for " 
           << aParticle->GetParticleName() << G4endl;
    return 0;
  }
  G4int nSubRegions = em->NumberOfSubCutoffRegions();

  if (1 < verbose) {
    G4cout << "G4LossTableManager::BuildTables() start to build range tables"
           << " and the sum of " << n_dedx << " processes"
           << " iem= " << iem << " em= " << em->GetProcessName()
           << " buildCSDARange= " << theParameters->BuildCSDARange()
           << " nSubRegions= " << nSubRegions;
    if(subcutProducer) { 
      G4cout << " SubCutProducer " << subcutProducer->GetName(); 
    }
    G4cout << G4endl;
  }
  // do not build tables if producer class is defined
  if(subcutProducer) { nSubRegions = 0; }

  dedx = em->DEDXTable(); 
  em->SetIonisation(true);
  em->SetDEDXTable(dedx, fIsIonisation);

  if (1 < n_dedx) {
    dedx = 0;
    dedx = G4PhysicsTableHelper::PreparePhysicsTable(dedx);
    tableBuilder->BuildDEDXTable(dedx, t_list);
    em->SetDEDXTable(dedx, fRestricted);
  }

  /*
  if(2==run && "e-" == aParticle->GetParticleName()) {
    G4cout << "G4LossTableManager::BuildTables for e- " << dedx << G4endl;
    G4cout << (*dedx) << G4endl;
    G4cout << "%%%%% Instance ID= " << (*dedx)[0]->GetInstanceID() << G4endl;
    G4cout << "%%%%% LastValue= " << (*dedx)[0]->GetLastValue() << G4endl;
    G4cout << "%%%%% 1.2 " << (*(dedx))[0]->Value(1.2) << G4endl;
  }
  */
  dedx_vector[iem] = dedx;

  G4PhysicsTable* range = em->RangeTableForLoss();
  if(!range) range  = G4PhysicsTableHelper::PreparePhysicsTable(range);
  range_vector[iem] = range;

  G4PhysicsTable* invrange = em->InverseRangeTable();
  if(!invrange) invrange = G4PhysicsTableHelper::PreparePhysicsTable(invrange);
  inv_range_vector[iem]  = invrange;

  tableBuilder->BuildRangeTable(dedx, range, true);
  tableBuilder->BuildInverseRangeTable(range, invrange, true);

  //  if(1<verbose) G4cout << *dedx << G4endl;

  em->SetRangeTableForLoss(range);
  em->SetInverseRangeTable(invrange);

  //  if(1<verbose) G4cout << *range << G4endl;

  std::vector<G4PhysicsTable*> listSub;
  std::vector<G4PhysicsTable*> listCSDA;

  for (i=0; i<n_dedx; ++i) {
    p = loss_list[i];
    if(p != em) { p->SetIonisation(false); }
    if(build_flags[i]) {
      p->SetLambdaTable(p->BuildLambdaTable(fRestricted));
    }
    if (0 < nSubRegions) {
      dedx = p->BuildDEDXTable(fSubRestricted);
      p->SetDEDXTable(dedx,fSubRestricted);
      listSub.push_back(dedx);
      if(build_flags[i]) {
        p->SetSubLambdaTable(p->BuildLambdaTable(fSubRestricted));
        if(p != em) { em->AddCollaborativeProcess(p); }
      }
    }
    if(theParameters->BuildCSDARange()) { 
      dedx = p->BuildDEDXTable(fTotal);
      p->SetDEDXTable(dedx,fTotal);
      listCSDA.push_back(dedx); 
    }     
  }

  if (0 < nSubRegions) {
    G4PhysicsTable* dedxSub = em->IonisationTableForSubsec();
    if (1 < listSub.size()) {
      em->SetDEDXTable(dedxSub, fIsSubIonisation);
      dedxSub = 0;
      dedxSub = G4PhysicsTableHelper::PreparePhysicsTable(dedxSub);
      tableBuilder->BuildDEDXTable(dedxSub, listSub);
      em->SetDEDXTable(dedxSub, fSubRestricted);
    }
  }
  if(theParameters->BuildCSDARange()) {
    G4PhysicsTable* dedxCSDA = em->DEDXunRestrictedTable();
    if (1 < n_dedx) {
      dedxCSDA = nullptr;
      dedxCSDA = G4PhysicsTableHelper::PreparePhysicsTable(dedxCSDA);
      tableBuilder->BuildDEDXTable(dedxCSDA, listCSDA);
      em->SetDEDXTable(dedxCSDA,fTotal);
    }
    G4PhysicsTable* rCSDA = em->CSDARangeTable();
    if(!rCSDA) { rCSDA = G4PhysicsTableHelper::PreparePhysicsTable(rCSDA); }
    tableBuilder->BuildRangeTable(dedxCSDA, rCSDA, true);
    em->SetCSDARangeTable(rCSDA);
  }

  if (1 < verbose) {
    G4cout << "G4LossTableManager::BuildTables: Tables are built for "
           << aParticle->GetParticleName()
           << "; ionisation process: " << em->GetProcessName()
           << "  " << em
           << G4endl;
  }
  return em;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::ParticleHaveNoLoss(
     const G4ParticleDefinition* aParticle)
{
  G4ExceptionDescription ed;
  ed << "Energy loss process not found for " << aParticle->GetParticleName() 
     << " !";
  G4Exception("G4LossTableManager::ParticleHaveNoLoss", "em0001",
              FatalException, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::SetVerbose(G4int val)
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VEnergyLossProcess*>& 
G4LossTableManager::GetEnergyLossProcessVector()
{
  return loss_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VEmProcess*>& G4LossTableManager::GetEmProcessVector()
{
  return emp_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

const std::vector<G4VMultipleScattering*>& 
G4LossTableManager::GetMultipleScatteringVector()
{
  return msc_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmSaturation* G4LossTableManager::EmSaturation()
{
  return theParameters->GetEmSaturation();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmConfigurator* G4LossTableManager::EmConfigurator()
{
  if(!emConfigurator) { emConfigurator = new G4EmConfigurator(verbose); }
  return emConfigurator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElectronIonPair* G4LossTableManager::ElectronIonPair()
{
  if(!emElectronIonPair) { 
    emElectronIonPair = new G4ElectronIonPair(verbose);
  }
  return emElectronIonPair;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
 
void G4LossTableManager::SetAtomDeexcitation(G4VAtomDeexcitation* p)
{
  if(atomDeexcitation != p) {
    delete atomDeexcitation;
    atomDeexcitation = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::SetSubCutProducer(G4VSubCutProducer* p) 
{
  if(subcutProducer != p) {
    delete subcutProducer;
    subcutProducer = p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableManager::PrintEWarning(G4String tit, G4double /*val*/)
{
  G4String ss = "G4LossTableManager::" + tit; 
  G4ExceptionDescription ed;
  /*
  ed << "Parameter is out of range: " << val 
     << " it will have no effect!\n" << " ## " 
     << " nbins= " << nbinsLambda 
     << " nbinsPerDecade= " << nbinsPerDecade 
     << " Emin(keV)= " << minKinEnergy/keV 
     << " Emax(GeV)= " << maxKinEnergy/GeV;
  */
  G4Exception(ss, "em0044", JustWarning, ed);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4LossTableManager::DumpHtml()
{
  // Automatic generation of html documentation page for physics lists
  // List processes and models for the most important
  // particles in descending order of importance

  char* dirName = getenv("G4PhysListDocDir");
  char* physListName = getenv("G4PhysListName");
  if (dirName && physListName) {
    G4String pathName = G4String(dirName) + "/" + G4String(physListName) 
                                          + ".html";

    std::ofstream outFile;
    outFile.open(pathName);
   
    outFile << "<html>\n";
    outFile << "<head>\n";
    outFile << "<title>Physics List Summary</title>\n";
    outFile << "</head>\n";
    outFile << "<body>\n";
    outFile << "<h2> Summary of Electromagnetic Processes, Models and Cross "
            << "Sections for Physics List "
            << G4String(physListName) << "</h2>\n";
    outFile << "<ul>\n";

    std::vector<G4ParticleDefinition*> particles {
        G4Gamma::Gamma(),
        G4Electron::Electron(),
        G4Positron::Positron(),
        G4Proton::ProtonDefinition(),
        G4MuonPlus::MuonPlusDefinition(),
        G4MuonMinus::MuonMinusDefinition(),
      };
   
    std::vector<G4VEmProcess*> emproc_vector = GetEmProcessVector();
    std::vector<G4VEnergyLossProcess*> enloss_vector = 
      GetEnergyLossProcessVector();
    std::vector<G4VMultipleScattering*> mscat_vector =
      GetMultipleScatteringVector();
    
    for (auto theParticle : particles) {
      outFile << "<li><h3>" << theParticle->GetParticleName() << "</h3><ul>\n";

			G4ProcessManager* pm = theParticle->GetProcessManager();
      G4ProcessVector*  pv = pm->GetProcessList();
      G4int plen = pm->GetProcessListLength();

      for (auto emproc : emproc_vector) {
        for (G4int i = 0; i < plen; ++i) {
          G4VProcess* proc = (*pv)[i];
          if (proc == emproc) {
            outFile << "<li>\n";
            proc->ProcessDescription(outFile);
            outFile << "</li>\n";
            break;
          }
        }
      }

      for (auto mscproc : mscat_vector) {
        for (G4int i = 0; i < plen; ++i) {
          G4VProcess* proc = (*pv)[i];
          if (proc == mscproc) {
            outFile << "<li>\n";
            proc->ProcessDescription(outFile);
            outFile << "</li>\n";
            break;
          }
        }
      }

      for (auto enlossproc : enloss_vector) {
        for (G4int i = 0; i < plen; ++i) {
          G4VProcess* proc = (*pv)[i];
          if (proc == enlossproc) {
            outFile << "<li>\n";
            proc->ProcessDescription(outFile);
            outFile << "</li>\n";
            break;
          }
        }
      }
      outFile << "</ul>\n";
      outFile << "</li>\n";
    }
    outFile << "</ul>\n";
    outFile << "</body>\n";
    outFile << "</html>\n";
    outFile.close();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

