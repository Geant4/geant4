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
// File name:     G4LossTableManager
//
// Author:        Vladimir Ivanchenko 
// 
// Creation date: 03.01.2002
//
// Modifications: 
//
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
//
// Class Description: 
//
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LossTableManager.hh"
#include "G4EnergyLossMessengerSTD.hh"
#include "G4PhysicsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4ProcessManager.hh"
#include "G4Electron.hh"

G4LossTableManager* G4LossTableManager::theInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager* G4LossTableManager::Instance()
{
  if(0 == theInstance) {
    theInstance = new G4LossTableManager();
  }
  return theInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::~G4LossTableManager()
{
  delete theMessenger;
  delete tableBuilder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4LossTableManager::G4LossTableManager()
{
  n_loss = 0;
  first_entry = true;
  all_tables_are_built = false;
  electron_table_are_built = false;
  currentLoss = 0;
  currentParticle = 0;
  lossFluctuationFlag = true;
  subCutoffFlag = false;
  rndmStepFlag = false;
  minSubRange = 0.0;
  maxRangeVariation = 1.0;
  maxFinalStep = 0.0;
  eIonisation = 0;
  minKinEnergy = 0.1*eV;
  maxKinEnergy = 100.0*GeV;
  theMessenger = new G4EnergyLossMessengerSTD();
  theElectron  = G4Electron::Electron();
  tableBuilder = new G4LossTableBuilder();
  integral = true;
  verbose = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

void G4LossTableManager::Clear()
{
  all_tables_are_built = false;
  electron_table_are_built = false;
  eIonisation = 0;
  currentLoss = 0;
  currentParticle = 0;
  if(n_loss) 
    {
      for(G4int i=0; i<n_loss; i++) 
        {
          if(dedx_vector[i]) dedx_vector[i]->clearAndDestroy();
          if(range_vector[i]) range_vector[i]->clearAndDestroy();
          if(inv_range_vector[i]) inv_range_vector[i]->clearAndDestroy();
	}
    
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

void G4LossTableManager::Register(G4VEnergyLossSTD* p)
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
  /*
  p->SetLossFluctuations(lossFluctuationFlag);
  p->SetSubCutoff(subCutoffFlag);
  p->SetRandomStep(rndmStepFlag);
  p->SetMinSubRange(minSubRange);
  p->SetMinKinEnergy(minKinEnergy);
  p->SetMaxKinEnergy(maxKinEnergy);
  p->SetStepLimits(maxRangeVariation, maxFinalStep);
  p->SetIntegral(integral);
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::Initialise(const G4ParticleDefinition* aParticle)
{
  if (first_entry) {
    first_entry = false;
    firstParticle = aParticle;
  }

  for(G4int j=0; j<n_loss; j++) {
    
    const G4ParticleDefinition* pd = loss_vector[j]->Particle();
    if (!pd) {
      pd = loss_vector[j]->GetProcessManager()->GetParticleType();
      loss_vector[j]->SetParticle(pd);
    }
    if(!tables_are_built[j]) all_tables_are_built = false;
    if(pd != part_vector[j]) {
      part_vector[j] = pd;
      base_part_vector[j] = loss_vector[j]->BaseParticle(); 
    }
    if(1 < verbose) {
      G4String nm = "unknown";
      if(pd) nm = pd->GetParticleName();
      G4cout << "For " << loss_vector[j]->GetProcessName()
             << " for " << nm
             << " tables_are_built= " << tables_are_built[j]
             << " procFlag= " << loss_vector[j]->TablesAreBuilt()
             << " all_tables_are_built= "     << all_tables_are_built  
             << G4endl;
    }
  } 

  // All tables have to be rebuilt
  if (all_tables_are_built && firstParticle == aParticle) {
    for (G4int i=0; i<n_loss; i++) {
      tables_are_built[i] = false;      
    }
    all_tables_are_built = false;
    electron_table_are_built = false;
    loss_map.clear();
  }

  if(1 < verbose) {
    G4cout << "electron_table_are_built= "  << electron_table_are_built
           << " all_tables_are_built= "     << all_tables_are_built  
           << G4endl;
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::BuildDEDXTable(const G4ParticleDefinition* aParticle)
{

  if(1 < verbose) {
    G4cout << "G4LossTableManager::BuildDEDXTable() for " 
           << aParticle->GetParticleName()
           << G4endl;
  }

  Initialise(aParticle);

  // identify all particles and built electron table
  if(!electron_table_are_built) {
    eIonisation = BuildTables(theElectron);
    electron_table_are_built = true;
  }

  // If particle is not an electron 
  if(aParticle != theElectron) {
    for(G4int i=0; i<n_loss; i++) {

      if(1 < verbose) {
        G4String nm = "0";
        G4String nm1= "0";
        G4String nm2= "0";
        const G4ParticleDefinition* pd = part_vector[i];
        const G4ParticleDefinition* bpd = base_part_vector[i];
        if(pd) nm = pd->GetParticleName();
        if(bpd) nm2 = bpd->GetParticleName();
        if(part_vector[i]) nm1 = part_vector[i]->GetParticleName();
        G4cout << "For " << loss_vector[i]->GetProcessName()
               << " for " << nm
	       << " (" << nm1 << ") "
               << " base_part= " << nm2 
               << " tables_are_built= " << tables_are_built[i]
               << G4endl;
      }
      if(aParticle == part_vector[i]) {
        if(tables_are_built[i]) break;
        if(!base_part_vector[i]) {
          G4VEnergyLossSTD* hIonisation = BuildTables(aParticle);
          for(G4int j=0; j<i; j++) {
            if(aParticle == base_part_vector[j]) {
              tables_are_built[j] = true;
              G4VEnergyLossSTD* em = loss_vector[j]; 
              em->Initialise();
              em->SetDEDXTable(hIonisation->DEDXTable());
              em->SetRangeTable(hIonisation->RangeTable());
              em->SetInverseRangeTable(hIonisation->InverseRangeTable());
              loss_map[part_vector[j]] = hIonisation;
              if(em->SecondaryParticle() == theElectron) {
                em->SetSecondaryRangeTable(eIonisation->RangeTable());
              }
            }
	  }
        } else {
          G4VEnergyLossSTD* em = loss_vector[i];
          for(G4int j=0; j<n_loss; j++) {
            if(part_vector[j] == base_part_vector[i]) {
              G4VEnergyLossSTD* hIonisation = loss_vector[j]; 
              if(tables_are_built[j] && hIonisation->DEDXTable()) {
                tables_are_built[i] = true;
                em->Initialise();
                em->SetDEDXTable(hIonisation->DEDXTable());
                em->SetRangeTable(hIonisation->RangeTable());
                em->SetInverseRangeTable(hIonisation->InverseRangeTable());
                loss_map[part_vector[i]] = hIonisation;
                if(em->SecondaryParticle() == theElectron) {
                  em->SetSecondaryRangeTable(eIonisation->RangeTable());
                }
                break;
              }
            }
          }
        }
        break;
      }
    }
  }

  all_tables_are_built = true;
  for (G4int i=0; i<n_loss; i++) {
    if ( !tables_are_built[i] ) {
      all_tables_are_built = false;
      break;
    }       
  }

  if(1 < verbose) {
    G4cout << "G4LossTableManager::BuildDEDXTable: end; " 
           << "all_tables_are_built= " << all_tables_are_built 
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::RetrieveDEDXTable(const G4ParticleDefinition* aParticle,
                                                 G4VEnergyLossSTD* theModel)
{
  for(G4int i=0; i<n_loss; i++) {
    if(theModel == loss_vector[i]) {
      tables_are_built[i] = true;
      if(theModel->DEDXTable()) {
        for(G4int j=0; j<n_loss; j++) {
          if(j != i && tables_are_built[j]) {
            if(aParticle == base_part_vector[j]) {
              G4VEnergyLossSTD* em = loss_vector[j]; 
              em->SetDEDXTable(theModel->DEDXTable());
              em->SetRangeTable(theModel->RangeTable());
              em->SetInverseRangeTable(theModel->InverseRangeTable());
              loss_map[part_vector[j]] = theModel;
	    }
	  }
	}
      } else {
        const G4ParticleDefinition* bp = theModel->BaseParticle();
        for(G4int j=0; j<n_loss; j++) {
          if(j != i && tables_are_built[j]) {
            if(bp == part_vector[j] && loss_vector[j]->DEDXTable()) {
              G4VEnergyLossSTD* hIonisation = loss_vector[j]; 
              theModel->SetDEDXTable(hIonisation->DEDXTable());
              theModel->SetRangeTable(hIonisation->RangeTable());
              theModel->SetInverseRangeTable(hIonisation->InverseRangeTable());
              loss_map[part_vector[j]] = hIonisation;
              break;
	    }
	  }
	}
      }
      break;
    }
  }

  // identify all particles and built electron table
  if(aParticle == theElectron) {
    electron_table_are_built = true;
    for(G4int i=0; i<n_loss; i++) {
      if(loss_vector[i]->SecondaryParticle() == theElectron) {
        loss_vector[i]->SetSecondaryRangeTable(theModel->RangeTable());
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VEnergyLossSTD* G4LossTableManager::BuildTables(const G4ParticleDefinition* aParticle)
{
  if(1 < verbose) {
    G4cout << "G4LossTableManager::BuildTables() for "
           << aParticle->GetParticleName() << G4endl;
  }

  // Check is it new particle or all tables have to be rebuilt
  G4std::vector<G4PhysicsTable*> list;
  list.clear();
  G4std::vector<G4VEnergyLossSTD*> loss_list;
  loss_list.clear();
  G4VEnergyLossSTD* em = 0; 
  G4int iem = 0;
  G4bool no_el = true;

  for(G4int i=0; i<n_loss; i++) {
    if(aParticle == part_vector[i]) {
      loss_vector[i]->Initialise();
      if(no_el && loss_vector[i]->SecondaryParticle() == theElectron) {
        em = loss_vector[i];
        iem= i;
        no_el = false;
      }
      if(!em) {
        em = loss_vector[i];
        iem= i;
      }
      tables_are_built[i] = true;
      list.push_back(loss_vector[i]->BuildDEDXTable());
      loss_list.push_back(loss_vector[i]);
    }
  }

  G4int n_dedx = list.size();
  if(!n_dedx) return 0;

  if(aParticle == theElectron) eIonisation = em;

  if(0 < verbose) {
    G4cout << "G4LossTableManager::BuildTables() start to build range tables"
           << " and the sum of " << n_dedx << " processes"
           << G4endl;
  }

  G4PhysicsTable* dedx = list[0];
  if(1 < n_dedx) {
    dedx = tableBuilder->BuildDEDXTable(list);
    for(G4int i=0; i<n_dedx; i++) {
      list[i]->clearAndDestroy();
    }
  }
  em->SetDEDXTable(dedx);
  dedx_vector[iem] = dedx;
  G4PhysicsTable* range = tableBuilder->BuildRangeTable(dedx);
  em->SetRangeTable(range);
  range_vector[iem] = range;
  G4PhysicsTable* invrange = tableBuilder->BuildInverseRangeTable(dedx, range);
  em->SetInverseRangeTable(invrange);
  inv_range_vector[iem] = invrange;

  loss_map[aParticle] = em;
  for(G4int j=0; j<n_dedx; j++) {
    if(loss_list[j]->SecondaryParticle() == theElectron) {
      loss_list[j]->SetSecondaryRangeTable(eIonisation->RangeTable());
    }
  }
  if(1 < verbose) {
    G4cout << "G4LossTableManager::BuildTables: Tables are built"
           << " for " << em->GetProcessName()
           << G4endl;
  }
  return em;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetLossFluctuations(G4bool val)
{
  lossFluctuationFlag = val;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetLossFluctuations(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetSubCutoff(G4bool val) 
{
  subCutoffFlag = val;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetSubCutoff(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetIntegral(G4bool val) 
{
  integral = val;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetIntegral(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMinSubRange(G4double val) 
{
  minSubRange = val;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetMinSubRange(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetRandomStep(G4bool val) 
{
  rndmStepFlag = val;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetRandomStep(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMinEnergy(G4double val) 
{
  minKinEnergy = val;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetMinKinEnergy(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetMaxEnergy(G4double val) 
{
  maxKinEnergy = val;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetMaxKinEnergy(val);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4LossTableManager::SetStepLimits(G4double v1, G4double v2) 
{
  maxRangeVariation = v1; 
  maxFinalStep = v2;
  for(G4int i=0; i<n_loss; i++) {
    loss_vector[i]->SetStepLimits(v1, v2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4LossTableManager::StorePhysicsTable(G4VEnergyLossSTD* el,
			   	       const G4String& dedx_file, 
			   	       const G4String& range_file, 
			   	       const G4String& inv_range_file, 
				             G4bool ascii)
{
  G4int i=0;
  for(i=0; i<n_loss; i++) {
    if(el == loss_vector[i]) break;
  }

  // store stopping power table
  if (dedx_vector[i]) {
    if ( !dedx_vector[i]->StorePhysicsTable(dedx_file, ascii) ){
      G4cout << "Fatal error theDEDXTable->StorePhysicsTable in <" 
             << dedx_file << ">"
             << G4endl;
      return false;
    }
  }

  if (range_vector[i]) {
    if ( !range_vector[i]->StorePhysicsTable(range_file, ascii) ){
      G4cout << "Fatal error theRangeTable->StorePhysicsTable in <" 
             << range_file << ">"
             << G4endl;
      return false;
    } 
  }

  if (inv_range_vector[i]) {
    if ( !inv_range_vector[i]->StorePhysicsTable(inv_range_file, ascii) ){
      G4cout << "Fatal error theInverseRangeTable->StorePhysicsTable in <" 
             << inv_range_file << ">"
             << G4endl;
      return false;
    }
  }
  return true;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


