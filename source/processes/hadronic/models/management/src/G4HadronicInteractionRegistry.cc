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
// $Id: G4HadronicInteractionRegistry.cc 97302 2016-06-01 09:30:11Z gcosmo $
//
// 23-Jan-2009 V.Ivanchenko make the class to be a singleton
// 17-Aug-2012 V.Ivanchenko added hadronic model factories

#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicInteraction.hh"

G4ThreadLocal G4HadronicInteractionRegistry* G4HadronicInteractionRegistry::instance = nullptr;

G4HadronicInteractionRegistry* G4HadronicInteractionRegistry::Instance()
{
  if(nullptr == instance) {
    static G4ThreadLocalSingleton<G4HadronicInteractionRegistry> inst;
    instance = inst.Instance();
  }
  return instance;
}

G4HadronicInteractionRegistry::G4HadronicInteractionRegistry()
  : isInitialized(false)
{
  //G4cout << "G4HadronicInteractionRegistry  " << this << G4endl;
}

G4HadronicInteractionRegistry::~G4HadronicInteractionRegistry()
{
  Clean();
}

void G4HadronicInteractionRegistry::Clean()
{
  size_t nModels = allModels.size();
  //G4cout << "G4HadronicInteractionRegistry::Clean() start " << nModels 
  //	    << " " << this << G4endl;
  for (size_t i=0; i<nModels; ++i) {
    if( allModels[i] ) {
      const char* xxx = (allModels[i]->GetModelName()).c_str();
      G4int len = (allModels[i]->GetModelName()).length();
      len = std::min(len, 9);
      const G4String mname = G4String(xxx, len);
      //G4cout << "G4HadronicInteractionRegistry: delete " << i << "  "
      //		<< allModels[i] << " " << mname 
      //		<< " " << this << G4endl;
      if(mname != "NeutronHP" && mname != "ParticleH") {
	delete allModels[i];
      }
      // G4cout << "done " << this << G4endl;
    }
  }
  allModels.clear();
  //G4cout <<"G4HadronicInteractionRegistry::Clean() is done "<<G4endl; 
}

void G4HadronicInteractionRegistry::InitialiseModels()
{
  size_t nModels = allModels.size();
  for (size_t i=0; i<nModels; ++i) {
    if( allModels[i] ) { allModels[i]->InitialiseModel(); }
  }
}

void 
G4HadronicInteractionRegistry::RegisterMe(G4HadronicInteraction * aModel)
{
  if(!aModel) { return; }
  size_t nModels = allModels.size();
  for (size_t i=0; i<nModels; ++i) {
    if( aModel == allModels[i] ) { return; }
  }
  //G4cout << "Register model <" << aModel->GetModelName() 
  //	 << ">  " << nModels+1 << "  " << aModel << G4endl;
  allModels.push_back(aModel);
}

void 
G4HadronicInteractionRegistry::RemoveMe(G4HadronicInteraction * aModel)
{
  if(!aModel) { return; }
  size_t nModels = allModels.size();
  for (size_t i=0; i<nModels; ++i) {
    if( aModel == allModels[i] ) {
      //G4cout << "DeRegister model <" << aModel 
      //	<< ">  " << i << G4endl;
      allModels[i] = 0;
      return;
    }
  }
}

G4HadronicInteraction* 
G4HadronicInteractionRegistry::FindModel(const G4String& name)
{
  G4HadronicInteraction* model = 0; 

  size_t nModels = allModels.size(); 
  for (size_t i=0; i<nModels; ++i) {
    G4HadronicInteraction* p = allModels[i]; 
    if(p) {
      if (p->GetModelName() == name) { 
	model = p;
	break; 
      }
    }
  }
  return model;
}

std::vector<G4HadronicInteraction*>
G4HadronicInteractionRegistry::FindAllModels(const G4String& name)
{
  std::vector<G4HadronicInteraction*> models;

  size_t nModels = allModels.size(); 
  for (size_t i=0; i<nModels; ++i) {
    G4HadronicInteraction* p = allModels[i]; 
    if(p) {
      if (p->GetModelName() == name) { 
        models.push_back(p);
      }
    }
  }
  return models;
}
