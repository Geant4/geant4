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
// $Id: G4EmConfigurator.cc,v 1.3 2008/11/21 12:30:29 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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

enum PType {unknown=0, eloss, discrete, msc};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmConfigurator::G4EmConfigurator()
{
  index = 1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
G4EmConfigurator::~G4EmConfigurator()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::AddExtraEmModel(const G4String& particleName,
				       G4VEmModel* em, 
				       G4VEmFluctuationModel* fm)
{
  particleList.push_back(particleName);
  modelList.push_back(em);
  flucModelList.push_back(fm);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::AddModelForRegion(const G4String& particleName,
					 const G4String& processName,
					 const G4String& modelName,
					 const G4String& regionName,
                                         G4double emin, G4double emax,
					 const G4String& flucModelName)
{
  particles.push_back(particleName);
  processes.push_back(processName);
  models.push_back(modelName);
  regions.push_back(regionName);
  flucModels.push_back(flucModelName);
  lowEnergy.push_back(emin);
  highEnergy.push_back(emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::SetExtraEmModel(const G4String& particleName,
				       const G4String& processName,
				       G4VEmModel* mod,
				       const G4String& regionName,
				       G4double emin,
				       G4double emax,
				       G4VEmFluctuationModel* fm)
{
  AddExtraEmModel(particleName, mod, fm);
  G4String fname = "";
  if(fm) fname = fm->GetName();
  AddModelForRegion(particleName, processName, mod->GetName(), regionName,
		    emin, emax, fname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::AddModels()
{
  size_t n = particles.size();
  //G4cout << " G4EmConfigurator::AddModels n= " << n << G4endl;
  if(n > 0) {
    for(size_t i=0; i<n; i++) {
      SetModelForRegion(particles[i],processes[i],models[i],regions[i],
                        flucModels[i],lowEnergy[i],highEnergy[i]);
    }
  }
  particles.clear();
  processes.clear();
  models.clear();
  flucModels.clear();
  regions.clear();
  lowEnergy.clear();
  highEnergy.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmConfigurator::SetModelForRegion(const G4String& particleName,
					 const G4String& processName,
					 const G4String& modelName,
					 const G4String& regionName,
					 const G4String& flucModelName,
					 G4double emin, G4double emax)
{
  //G4cout << " G4EmConfigurator::SetModelForRegion" << G4endl;

  // new set
  index--;

  G4ParticleTable::G4PTblDicIterator* theParticleIterator = 
    G4ParticleTable::GetParticleTable()->GetIterator(); 

  theParticleIterator->reset();
  while( (*theParticleIterator)() ) {
    const G4ParticleDefinition* part = theParticleIterator->value();

    //G4cout << particleName << " " << part->GetParticleName() << G4endl;

    if(particleName == part->GetParticleName() || 
       (particleName == "charged" && part->GetPDGCharge() != 0.0) ) {


      // search for process
      G4ProcessManager* pmanager = part->GetProcessManager();
      G4ProcessVector* plist = pmanager->GetProcessList();
      G4int np = pmanager->GetProcessListLength();
  
      //G4cout << processName << " in list of " << np << G4endl;

      G4VProcess* proc = 0;
      for(G4int i=0; i<np; i++) {
	if(processName == (*plist)[i]->GetProcessName()) {
	  proc = (*plist)[i];
	  break;
	}
      }
      if(!proc) {
	G4cout << "### G4EmConfigurator WARNING: fails to find a process <"
	       << processName << "> for " << particleName << G4endl;
	
      } else {

	// classify process
	PType ptype = discrete;
	G4int ii = proc->GetProcessSubType();
	if(10 == ii) ptype = msc;
	else if(2 <= ii && 4 >= ii) ptype = eloss;

	// find out model     
	G4VEmModel* mod = 0;
	G4VEmFluctuationModel* fluc = 0;

	G4int nm = modelList.size();
	//G4cout << "Search model " << modelName << " in " << nm << G4endl;

	for(G4int i=0; i<nm; i++) {
	  if(modelName == modelList[i]->GetName() &&
             (particleList[i] == "" || particleList[i] == particleName) ) {
	    mod  = modelList[i];
	    fluc = flucModelList[i];
	    break;
	  }
	}

	if("dummy" == modelName) mod = new G4DummyModel();

	if(!mod) {
	  G4cout << "### G4EmConfigurator WARNING: fails to find a model <"
		 << modelName << "> for process <" 
		 << processName << "> and " << particleName 
		 << G4endl;
	  if(flucModelName != "")  
	    G4cout << "                            fluctuation model <" 
		   << flucModelName << G4endl;
	} else {

	  // search for region
	  G4Region* reg = 0;
	  G4RegionStore* regStore = G4RegionStore::GetInstance();
	  G4String r = regionName;
	  if(r == "" || r == "world" || r == "World") r = "DefaultRegionForTheWorld";
	  reg = regStore->GetRegion(r, true); 
	  if(!reg) {
	    G4cout << "### G4EmConfigurator WARNING: fails to find a region <"
		   << r << "> for model <" << modelName << "> of the process " 
		   << processName << " and " << particleName << G4endl;
	    return;
	  }

	  // energy limits
	  G4double e1 = std::max(emin,mod->LowEnergyLimit());
	  G4double e2 = std::min(emax,mod->HighEnergyLimit());
	  if(e2 < e1) e2 = e1;
	  mod->SetLowEnergyLimit(e1);
	  mod->SetHighEnergyLimit(e2);

	  //G4cout << "e1= " << e1 << " e2= " << e2 << G4endl;

	  // added model
	  if(ptype == eloss) {
	    G4VEnergyLossProcess* p = reinterpret_cast<G4VEnergyLossProcess*>(proc);
	    p->AddEmModel(index,mod,fluc,reg);
	    //G4cout << "### Added eloss model order= " << index << " for " 
	    //	   << particleName << " and " << processName << "  " << mod << G4endl;
	  } else if(ptype == discrete) {
	    G4VEmProcess* p = reinterpret_cast<G4VEmProcess*>(proc);
	    p->AddEmModel(index,mod,reg);
	  } else if(ptype == msc) {
	    G4VMultipleScattering* p = reinterpret_cast<G4VMultipleScattering*>(proc);
	    p->AddEmModel(index,mod,reg);
	  }
	}
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......








