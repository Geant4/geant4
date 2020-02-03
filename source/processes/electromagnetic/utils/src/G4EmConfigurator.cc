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
// -------------------------------------------------------------------
//
// GEANT4 Class 
//
// File name:     G4EmConfigurator
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 14.07.2008
//
// Modifications:
//
// Class Description:
//
// This class provides configuration EM models for 
// particles/processes/regions
//

// -------------------------------------------------------------------
//

#include "G4EmConfigurator.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4VProcess.hh"
#include "G4ProcessVector.hh"
#include "G4RegionStore.hh"
#include "G4Region.hh"
#include "G4DummyModel.hh"
#include "G4VEnergyLossProcess.hh"
#include "G4VEmProcess.hh"
#include "G4VMultipleScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmConfigurator::G4EmConfigurator(G4int val):verbose(val)
{
  index = -10;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4EmConfigurator::~G4EmConfigurator()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::SetExtraEmModel(const G4String& particleName,
                                       const G4String& processName,
                                       G4VEmModel* mod,
                                       const G4String& regionName,
                                       G4double emin,
                                       G4double emax,
                                       G4VEmFluctuationModel* fm)
{
  if(!mod) { return; }
  if(1 < verbose) {
    G4cout << " G4EmConfigurator::SetExtraEmModel " << mod->GetName()
           << " for " << particleName 
           << " and " << processName 
           << " in the region <" << regionName
           << "> Emin(MeV)= " << emin/MeV
           << " Emax(MeV)= " << emax/MeV
           << G4endl;
  }

  models.push_back(mod);
  flucModels.push_back(fm);
  G4double emin0 = std::max(emin, mod->LowEnergyLimit());
  G4double emax0 = std::min(emax, mod->HighEnergyLimit());
  mod->SetActivationHighEnergyLimit(emax0);

  particles.push_back(particleName);
  processes.push_back(processName);
  regions.push_back(regionName);
  lowEnergy.push_back(emin0);
  highEnergy.push_back(emax0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::AddModels()
{
  size_t n = models.size();
  if(0 < verbose) {
    G4cout << "### G4EmConfigurator::AddModels n= " << n << G4endl;
  }
  if(n > 0) {
    for(size_t i=0; i<n; ++i) {
      if(models[i]) {
        G4Region* reg = FindRegion(regions[i]);
        if(reg) {
          --index;
          SetModelForRegion(models[i],flucModels[i],reg,
                            particles[i],processes[i],
                            lowEnergy[i],highEnergy[i]);
        }
      }
    }
  }
  Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::SetModelForRegion(G4VEmModel* mod,
                                         G4VEmFluctuationModel* fm,
                                         G4Region* reg,
                                         const G4String& particleName,
                                         const G4String& processName,
                                         G4double emin, G4double emax)
{
  if(!mod) { return; }
  if(1 < verbose) {
    G4cout << " G4EmConfigurator::SetModelForRegion: " << mod->GetName() 
           << G4endl;
    G4cout << " For " << particleName 
           << " and " << processName 
           << " in the region <" << reg->GetName()
           << " Emin(MeV)= " << emin/MeV
           << " Emax(MeV)= " << emax/MeV;
    if(fm) { G4cout << " FLmodel " << fm->GetName(); }
    G4cout << G4endl;
  }

  // Loop checking, 03-Aug-2015, Vladimir Ivanchenko
  auto myParticleIterator = G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() ) {
    const G4ParticleDefinition* part = myParticleIterator->value();

    if((part->GetParticleName() == particleName) ||
       (particleName == "all") ||
       (particleName == "charged" && part->GetPDGCharge() != 0.0)) {

      // search for process
      G4ProcessManager* pmanager = part->GetProcessManager();
      G4ProcessVector* plist = pmanager->GetProcessList();
      G4int np = pmanager->GetProcessListLength();
  
      if(1 < verbose) {
	G4cout << "Check process <" << processName << "> for " 
	       << particleName << " in list of " << np << " processes" 
	       << G4endl;
      }
      G4VProcess* proc = nullptr;
      for(G4int i=0; i<np; ++i) {
        if(processName == (*plist)[i]->GetProcessName()) {
          proc = (*plist)[i];
          break;
        }
      }
      if(!proc) {
        G4cout << "### G4EmConfigurator WARNING: fails to find a process <"
               << processName << "> for " << particleName << G4endl;
        return;        
      } 

      if(!UpdateModelEnergyRange(mod, emin, emax)) { return; }
      // classify process
      G4int ii = proc->GetProcessSubType();
      if(10 == ii) {
        G4VMultipleScattering* p = static_cast<G4VMultipleScattering*>(proc);
	p->AddEmModel(index,mod,reg);
	if(1 < verbose) {
	  G4cout << "### Added msc model order= " << index << " for " 
		 << particleName << " and " << processName << G4endl;
	}
      } else if(2 <= ii && 4 >= ii) {
        G4VEnergyLossProcess* p = static_cast<G4VEnergyLossProcess*>(proc);
	p->AddEmModel(index,mod,fm,reg);
	if(1 < verbose) {
	  G4cout << "### Added eloss model order= " << index << " for " 
		 << particleName << " and " << processName << G4endl;
        }
      } else {
        G4VEmProcess* p = static_cast<G4VEmProcess*>(proc);
        p->AddEmModel(index,mod,reg);
        if(1 < verbose) {
          G4cout << "### Added em model order= " << index << " for " 
                 << particleName << " and " << processName << G4endl;
        }
      } 
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmConfigurator::PrepareModels(const G4ParticleDefinition* aParticle,
                                G4VEnergyLossProcess* p)
{
  size_t n = particles.size();
  if(1 < verbose) {
    G4cout << " G4EmConfigurator::PrepareModels for EnergyLoss n= " 
           << n << G4endl;
  }
  if(n > 0) {
    G4String particleName = aParticle->GetParticleName(); 
    G4String processName  = p->GetProcessName(); 
    //G4cout <<  particleName << "  " <<  processName << G4endl;
    for(size_t i=0; i<n; ++i) {
      //G4cout <<  particles[i] << "  " <<  processes[i] << G4endl;
      if(processName == processes[i]) {
        if((particleName == particles[i]) ||
           (particles[i] == "all") ||
           (particles[i] == "charged" && aParticle->GetPDGCharge() != 0.0)) {
          G4Region* reg = FindRegion(regions[i]);
          //G4cout << "Region " << reg << G4endl;
          if(reg) {
            --index;
            G4VEmModel* mod = models[i];
            G4VEmFluctuationModel* fm = flucModels[i];
            if(mod) {
              if(UpdateModelEnergyRange(mod, lowEnergy[i], highEnergy[i])) {
                p->AddEmModel(index,mod,fm,reg);
                if(1 < verbose) {
                  G4cout << "### Added eloss model order= " << index << " for " 
                         << particleName << " and " << processName << G4endl;
                }
              }
            } else if(fm) {
              p->SetFluctModel(fm);
            }
          }
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmConfigurator::PrepareModels(const G4ParticleDefinition* aParticle,
                                G4VEmProcess* p)
{
  size_t n = particles.size();
  if(1 < verbose) {
    G4cout << " G4EmConfigurator::PrepareModels for EM process n= " 
           << n << G4endl;
  }
  if(n > 0) {
    G4String particleName = aParticle->GetParticleName(); 
    G4String processName  = p->GetProcessName(); 
    //G4cout <<  particleName << "  " <<  particleName << G4endl;
    for(size_t i=0; i<n; ++i) {
      if(processName == processes[i]) {
        if((particleName == particles[i]) ||
           (particles[i] == "all") ||
           (particles[i] == "charged" && aParticle->GetPDGCharge() != 0.0)) {
          G4Region* reg = FindRegion(regions[i]);
          //G4cout << "Region " << reg << G4endl;
          if(reg) {
            --index;
            G4VEmModel* mod = models[i];
            if(mod) {
              if(UpdateModelEnergyRange(mod, lowEnergy[i], highEnergy[i])) { 
                p->AddEmModel(index,mod,reg);
                if(1 < verbose) {
                  G4cout << "### Added em model order= " << index << " for " 
                         << particleName << " and " << processName << G4endl;
                }
              }
            }
          }
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void 
G4EmConfigurator::PrepareModels(const G4ParticleDefinition* aParticle,
                                G4VMultipleScattering* p)
{
  size_t n = particles.size();
  if(1 < verbose) {
    G4cout << " G4EmConfigurator::PrepareModels for MSC process n= " 
           << n << G4endl;
  }

  if(n > 0) {
    G4String particleName = aParticle->GetParticleName(); 
    G4String processName  = p->GetProcessName(); 
    for(size_t i=0; i<n; ++i) {
      if(processName == processes[i]) {
        if((particleName == particles[i]) ||
           (particles[i] == "all") ||
           (particles[i] == "charged" && aParticle->GetPDGCharge() != 0.0)) {
          G4Region* reg = FindRegion(regions[i]);
          if(reg) {
            --index;
            G4VEmModel* mod = models[i];
            if(mod) {
              if(UpdateModelEnergyRange(mod, lowEnergy[i], highEnergy[i])) { 
                p->AddEmModel(index,mod,reg);
                //G4cout << "### Added msc model order= " << index << " for " 
                //     << particleName << " and " << processName << G4endl;
              }
            }
          }
        }
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::Clear()
{
  particles.clear();
  processes.clear();
  models.clear();
  flucModels.clear();
  regions.clear();
  lowEnergy.clear();
  highEnergy.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Region* G4EmConfigurator::FindRegion(const G4String& regionName)
{
  // search for region
  G4Region* reg = 0;
  G4RegionStore* regStore = G4RegionStore::GetInstance();
  G4String r = regionName;
  if(r == "" || r == "world" || r == "World") { 
    r = "DefaultRegionForTheWorld";
  }
  reg = regStore->GetRegion(r, true); 
  if(!reg) {
    G4cout << "### G4EmConfigurator WARNING: fails to find a region <"
           << r << G4endl;
  } else if(verbose > 1) {
    G4cout << "### G4EmConfigurator finds out G4Region <" << r << ">" 
           << G4endl;
  } 
  return reg;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmConfigurator::UpdateModelEnergyRange(G4VEmModel* mod, 
                                                G4double emin, G4double emax)
{
  // energy limits
  G4double e1 = std::max(emin,mod->LowEnergyLimit());
  G4double e2 = std::min(emax,mod->HighEnergyLimit());
  if(e2 <= e1) {
    G4cout << "### G4EmConfigurator WARNING: empty energy interval"
           << " for <" << mod->GetName() 
           << ">  Emin(MeV)= " << e1/CLHEP::MeV 
           << ">  Emax(MeV)= " << e2/CLHEP::MeV 
           << G4endl;
    return false;        
  }
  mod->SetLowEnergyLimit(e1);
  mod->SetHighEnergyLimit(e2);
  if(verbose > 1) {
    G4cout << "### G4EmConfigurator for " << mod->GetName() 
           << " Emin(MeV)= " << e1/MeV << " Emax(MeV)= " << e2/MeV 
           << G4endl;
  } 
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
