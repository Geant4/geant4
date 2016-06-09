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
// $Id: G4LossTableManager.cc,v 1.51 2004/12/09 10:38:02 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
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
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
//
// Class Description:
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LossTableManager.hh"
#include "G4EnergyLossMessenger.hh"
#include "G4PhysicsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4VMultipleScattering.hh"
#include "G4VEmProcess.hh"
#include "G4ProductionCutsTable.hh"
#include "G4PhysicsTableHelper.hh"

G4LossTableManager* G4LossTableManager::theInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager* G4LossTableManager::Instance()
{
  if(0 == theInstance) {
    static G4LossTableManager manager;
    theInstance = &manager;
  }
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::~G4LossTableManager()
{
  for (G4int i=0; i<n_loss; i++) {
    if( loss_vector[i] ) delete loss_vector[i];
  }
  size_t msc = msc_vector.size();
  for (size_t j=0; j<msc; j++) {
    if(msc_vector[j] ) delete msc_vector[j];
  }
  size_t emp = emp_vector.size();
  for (size_t k=0; k<emp; k++) {
    if(emp_vector[k] ) delete emp_vector[k];
  }
  Clear();
  delete theMessenger;
  delete tableBuilder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::G4LossTableManager()
{
  n_loss = 0;
  first_entry = true;
  all_tables_are_built = false;
  all_tables_are_stored = false;
  currentLoss = 0;
  currentParticle = 0;
  lossFluctuationFlag = true;
  subCutoffFlag = false;
  rndmStepFlag = false;
  minSubRange = 0.0;
  maxRangeVariation = 1.0;
  maxFinalStep = 0.0;
  minKinEnergy = 0.1*keV;
  maxKinEnergy = 100.0*TeV;
  maxKinEnergyForMuons = 100.*TeV;
  theMessenger = new G4EnergyLossMessenger();
  theElectron  = G4Electron::Electron();
  tableBuilder = new G4LossTableBuilder();
  integral = true;
  integralActive = false;
  buildPreciseRange = false;
  minEnergyActive = false;
  maxEnergyActive = false;
  maxEnergyForMuonsActive = false;
  stepFunctionActive = false;
  verbose = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Clear()
{
  all_tables_are_built = false;
  currentLoss = 0;
  currentParticle = 0;
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
      n_loss = 0;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEnergyLossProcess* p)
{
  n_loss++;
  loss_vector.push_back(p);
  part_vector.push_back(0);
  base_part_vector.push_back(0);
  dedx_vector.push_back(0);
  range_vector.push_back(0);
  inv_range_vector.push_back(0);
  tables_are_built.push_back(false);
  all_tables_are_built = false;
  if(!lossFluctuationFlag) p->SetLossFluctuations(false);
  if(subCutoffFlag)        p->SetSubCutoff(true);
  if(stepFunctionActive)   p->SetStepFunction(maxRangeVariation, maxFinalStep);
  if(integralActive)       p->SetIntegral(integral);
  if(minEnergyActive)      p->SetMinKinEnergy(minKinEnergy);
  if(maxEnergyActive)      p->SetMaxKinEnergy(maxKinEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEnergyLossProcess* p)
{
  for (G4int i=0; i<n_loss; i++) {
    if(loss_vector[i] == p) loss_vector[i] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VMultipleScattering* p)
{
  msc_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VMultipleScattering* p)
{
  size_t msc = msc_vector.size();
  for (size_t i=0; i<msc; i++) {
    if(msc_vector[i] == p) msc_vector[i] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Register(G4VEmProcess* p)
{
  emp_vector.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::DeRegister(G4VEmProcess* p)
{
  size_t emp = emp_vector.size();
  for (size_t i=0; i<emp; i++) {
    if(emp_vector[i] == p) emp_vector[i] = 0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::RegisterIon(const G4ParticleDefinition* ion, G4VEnergyLossProcess* p)
{
  loss_map[ion] = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::RegisterExtraParticle(const G4ParticleDefinition* part,
                                                     G4VEnergyLossProcess* p)
{
  n_loss++;
  loss_vector.push_back(p);
  part_vector.push_back(part);
  base_part_vector.push_back(p->BaseParticle());
  dedx_vector.push_back(0);
  range_vector.push_back(0);
  inv_range_vector.push_back(0);
  tables_are_built.push_back(false);
  all_tables_are_built = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::EnergyLossProcessIsInitialised(
                   const G4ParticleDefinition* particle, G4VEnergyLossProcess* p)
{
  if (first_entry || (particle == firstParticle && all_tables_are_built) ) {
    all_tables_are_built = true;
    loss_map.clear();
    for (G4int i=0; i<n_loss; i++) {
      G4VEnergyLossProcess* el = loss_vector[i];

      if(el) {

	const G4ProcessManager* pm = el->GetProcessManager();
        if (pm->GetProcessActivation(el)) tables_are_built[i] = false;
        else                              tables_are_built[i] = true;

	if (!tables_are_built[i]) all_tables_are_built = false;

      } else {
        tables_are_built[i] = true;
        part_vector[i] = 0;

      }

    }

    if (first_entry) {
      first_entry = false;
      firstParticle = particle;
    }
  }

  SetParameters(p);
  for (G4int j=0; j<n_loss; j++) {
    if (p == loss_vector[j]) {

      if (!part_vector[j]) {
        part_vector[j] = particle;
      	base_part_vector[j] = p->BaseParticle();
      }
      if(maxEnergyForMuonsActive) {
        G4double dm = std::fabs(particle->GetPDGMass() - 105.7*MeV);
	if(dm < 5.*MeV) p->SetMaxKinEnergy(maxKinEnergyForMuons);
      }

      if(0 < verbose) {
        G4cout << "For " << p->GetProcessName()
               << " for " << part_vector[j]->GetParticleName()
               << " tables_are_built= " << tables_are_built[j]
               << " procFlag= " << loss_vector[j]->TablesAreBuilt()
               << " all_tables_are_built= "     << all_tables_are_built
               << G4endl;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EnergyLossMessenger* G4LossTableManager::GetMessenger()
{
  return theMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::ParticleHaveNoLoss(const G4ParticleDefinition* aParticle)
{
  G4String s = "G4LossTableManager:: dE/dx table not found for "
             + aParticle->GetParticleName() + "!";
  G4Exception(s);
  exit(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4LossTableManager::BuildPreciseRange() const
{
  return buildPreciseRange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::BuildPhysicsTable(const G4ParticleDefinition* aParticle,
                                                 G4VEnergyLossProcess* p)
{

  if(0 < verbose) {
    G4cout << "### G4LossTableManager::BuildDEDXTable() is requested for "
           << aParticle->GetParticleName()
	   << " and process " << p->GetProcessName()
           << G4endl;
  }
  if (all_tables_are_built) return;
  all_tables_are_built = true;

  for(G4int i=0; i<n_loss; i++) {

    if(!tables_are_built[i] && !base_part_vector[i]) {

      const G4ParticleDefinition* curr_part = part_vector[i];
      G4VEnergyLossProcess* curr_proc = BuildTables(curr_part);
      CopyTables(curr_part, curr_proc);
    }
  }

  for (G4int ii=0; ii<n_loss; ii++) {
    if ( !tables_are_built[ii] ) {
      all_tables_are_built = false;
      break;
    }
  }

  if(0 < verbose) {
    G4cout << "### G4LossTableManager::BuildDEDXTable end: "
           << "all_tables_are_built= " << all_tables_are_built
           << G4endl;

    if(all_tables_are_built)
      G4cout << "### All dEdx and Range tables are built #####" << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::CopyTables(const G4ParticleDefinition* part,
                                          G4VEnergyLossProcess* base_proc)
{
  for (G4int j=0; j<n_loss; j++) {

    G4VEnergyLossProcess* proc = loss_vector[j];
    if(proc == base_proc || proc->Particle() == part) tables_are_built[j] = true;

    if (!tables_are_built[j] && part == base_part_vector[j]) {
      tables_are_built[j] = true;
      proc->SetDEDXTable(base_proc->DEDXTable());
      proc->SetDEDXunRestrictedTable(base_proc->DEDXunRestrictedTable());
      proc->SetPreciseRangeTable(base_proc->PreciseRangeTable());
      proc->SetRangeTableForLoss(base_proc->RangeTableForLoss());
      proc->SetInverseRangeTable(base_proc->InverseRangeTable());
      proc->SetLambdaTable(base_proc->LambdaTable());
      proc->SetSubLambdaTable(base_proc->SubLambdaTable());
      loss_map[part_vector[j]] = proc;
      if (0 < verbose) {
         G4cout << "For " << proc->GetProcessName()
                << " for " << part_vector[j]->GetParticleName()
                << " base_part= " << part->GetParticleName()
                << " tables are assigned "
                << G4endl;
      }
    }

    if (theElectron == part && theElectron == proc->SecondaryParticle() )
      proc->SetSecondaryRangeTable(base_proc->RangeTableForLoss());
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEnergyLossProcess* G4LossTableManager::BuildTables(const G4ParticleDefinition* aParticle)
{
  if(0 < verbose) {
    G4cout << "G4LossTableManager::BuildTables() for "
           << aParticle->GetParticleName() << G4endl;
  }

  std::vector<G4PhysicsTable*> list;
  list.clear();
  std::vector<G4VEnergyLossProcess*> loss_list;
  loss_list.clear();
  G4VEnergyLossProcess* em = 0;
  G4int iem = 0;

  for (G4int i=0; i<n_loss; i++) {
    if (aParticle == part_vector[i] && !tables_are_built[i] && loss_vector[i]) {
      if (loss_vector[i]->IsIonisationProcess()) {
        em = loss_vector[i];
        iem= i;
      }
      list.push_back(loss_vector[i]->BuildDEDXTable());
      loss_list.push_back(loss_vector[i]);
      tables_are_built[i] = true;
    }
  }

  G4int n_dedx = list.size();
  if (!n_dedx) return 0;

  if (1 < verbose) {
    G4cout << "G4LossTableManager::BuildTables() start to build range tables"
           << " and the sum of " << n_dedx << " processes"
           << " iem= " << iem << " em= " << em
	   << G4endl;
  }

  G4PhysicsTable* dedx = em->DEDXTable();
  if (1 < n_dedx) tableBuilder->BuildDEDXTable(dedx, list);

  dedx_vector[iem] = dedx;

  G4PhysicsTable* range = em->RangeTableForLoss();
  if(!range) range  = G4PhysicsTableHelper::PreparePhysicsTable(range);
  range_vector[iem] = range;

  G4PhysicsTable* invrange = em->InverseRangeTable();
  if(!invrange) invrange = G4PhysicsTableHelper::PreparePhysicsTable(invrange);
  inv_range_vector[iem]  = invrange;

  tableBuilder->BuildRangeTable(dedx, range);
  tableBuilder->BuildInverseRangeTable(range, invrange);

  if(buildPreciseRange) {
    std::vector<G4PhysicsTable*> newlist;
    newlist.clear();
    G4PhysicsTable* dedxForRange = em->DEDXunRestrictedTable();
    for (G4int i=0; i<n_dedx; i++) {
      newlist.push_back(loss_list[i]->BuildDEDXTableForPreciseRange());
    }
    if (1 < n_dedx) tableBuilder->BuildDEDXTable(dedxForRange, newlist);

    em->SetDEDXunRestrictedTable(dedxForRange);
    range = em->PreciseRangeTable();
    if(!range) range  = G4PhysicsTableHelper::PreparePhysicsTable(range);
    tableBuilder->BuildRangeTable(dedxForRange, range);
    em->SetPreciseRangeTable(range);
  }

  loss_map[aParticle] = em;
  for (G4int j=0; j<n_dedx; j++) {
    G4PhysicsTable* lambdaTable = loss_list[j]->BuildLambdaTable();
    loss_list[j]->SetLambdaTable(lambdaTable);
    if (0 < loss_list[j]->NumberOfSubCutoffRegions()) {
      lambdaTable = loss_list[j]->BuildLambdaSubTable();
      loss_list[j]->SetSubLambdaTable(lambdaTable);
    }
  }
  if (0 < verbose) {
    G4cout << "G4LossTableManager::BuildTables: Tables are built for "
           << aParticle->GetParticleName()
	   << "; ionisation process: " << em->GetProcessName()
           << G4endl;
  }
  return em;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetLossFluctuations(G4bool val)
{
  lossFluctuationFlag = val;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetLossFluctuations(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetSubCutoff(G4bool val)
{
  subCutoffFlag = val;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetSubCutoff(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetIntegral(G4bool val)
{
  integral = val;
  integralActive = true;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetIntegral(val);
  }
  size_t emp = emp_vector.size();
  for (size_t k=0; k<emp; k++) {
    if(emp_vector[k]) emp_vector[k]->SetIntegral(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMinSubRange(G4double val)
{
  minSubRange = val;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetMinSubRange(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetRandomStep(G4bool)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMinEnergy(G4double val)
{
  minEnergyActive = true;
  minKinEnergy = val;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetMinKinEnergy(val);
  }
  size_t msc = msc_vector.size();
  for (size_t j=0; j<msc; j++) {
    if(msc_vector[j]) msc_vector[j]->SetMinKinEnergy(val);
  }
  size_t emp = emp_vector.size();
  for (size_t k=0; k<emp; k++) {
    if(emp_vector[k]) emp_vector[k]->SetMinKinEnergy(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMaxEnergy(G4double val)
{
  maxEnergyActive = true;
  maxKinEnergy = val;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetMaxKinEnergy(val);
  }
  size_t msc = msc_vector.size();
  for (size_t j=0; j<msc; j++) {
    if(msc_vector[j]) msc_vector[j]->SetMaxKinEnergy(val);
  }
  size_t emp = emp_vector.size();
  for (size_t k=0; k<emp; k++) {
    if(emp_vector[k]) emp_vector[k]->SetMaxKinEnergy(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMaxEnergyForPreciseRange(G4double val)
{
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetMaxKinEnergyForPreciseRange(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMaxEnergyForMuons(G4double val)
{
  maxEnergyForMuonsActive = true;
  maxKinEnergyForMuons = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetDEDXBinning(G4int val)
{
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetDEDXBinning(val);
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetDEDXBinningForPreciseRange(G4int val)
{
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetDEDXBinningForPreciseRange(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetLambdaBinning(G4int val)
{
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetLambdaBinning(val);
  }
  size_t msc = msc_vector.size();
  for (size_t j=0; j<msc; j++) {
    if(msc_vector[j]) msc_vector[j]->SetBinning(val);
  }
  size_t emp = emp_vector.size();
  for (size_t k=0; k<emp; k++) {
    if(emp_vector[k]) emp_vector[k]->SetLambdaBinning(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetVerbose(G4int val)
{
  verbose = val;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetVerboseLevel(val);
  }
  size_t msc = msc_vector.size();
  for (size_t j=0; j<msc; j++) {
    if(msc_vector[j]) msc_vector[j]->SetVerboseLevel(val);
  }
  size_t emp = emp_vector.size();
  for (size_t k=0; k<emp; k++) {
    if(emp_vector[k]) emp_vector[k]->SetVerboseLevel(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetStepLimits(G4double v1, G4double v2)
{
  stepFunctionActive = true;
  maxRangeVariation = v1;
  maxFinalStep = v2;
  for(G4int i=0; i<n_loss; i++) {
    if(loss_vector[i]) loss_vector[i]->SetStepFunction(v1, v2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetBuildPreciseRange(G4bool val)
{
  buildPreciseRange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetParameters(G4VEnergyLossProcess* p)
{
  if(stepFunctionActive) p->SetStepFunction(maxRangeVariation, maxFinalStep);
  if(integralActive)     p->SetIntegral(integral);
  if(minEnergyActive)    p->SetMinKinEnergy(minKinEnergy);
  if(maxEnergyActive)    p->SetMaxKinEnergy(maxKinEnergy);
  p->SetVerboseLevel(verbose);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4VEnergyLossProcess*>& G4LossTableManager::GetEnergyLossProcessVector()
{
  return loss_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4VEmProcess*>& G4LossTableManager::GetEmProcessVector()
{
  return emp_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::vector<G4VMultipleScattering*>& G4LossTableManager::GetMultipleScatteringVector()
{
  return msc_vector;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
