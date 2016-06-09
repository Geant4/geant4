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
// $Id: G4HadronicInteractionRegistry.cc,v 1.3.2.1 2009/03/03 11:26:45 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-02 $
//
// 23-Jan-2009 V.Ivanchenko make the class to be a singleton

#include "G4HadronicInteractionRegistry.hh"
#include "G4HadronicInteraction.hh"

G4HadronicInteractionRegistry* G4HadronicInteractionRegistry::theInstance = 0;

G4HadronicInteractionRegistry* G4HadronicInteractionRegistry::Instance()
{
  if(0 == theInstance) {
    static G4HadronicInteractionRegistry manager;
    theInstance = &manager;
  }
  return theInstance;
}

G4HadronicInteractionRegistry::G4HadronicInteractionRegistry()
{
  nModels = 0;
}

G4HadronicInteractionRegistry::~G4HadronicInteractionRegistry()
{
  Clean();
}

void G4HadronicInteractionRegistry::Clean()
{
  //G4cout << "G4HadronicInteractionRegistry::Clean() start " << nModels << G4endl;
  if(0 < nModels) {
    for (G4int i=0; i<nModels; i++) {
      if( allModels[i] ) {
	//G4cout << "delete " << i << G4endl;
        //G4cout << allModels[i]->GetModelName() << G4endl;
	delete allModels[i];
	allModels[i] = 0;
      }
    }
  }
  //G4cout << "G4HadronicInteractionRegistry::Clean() is done " << G4endl; 
  nModels = 0;
}

void G4HadronicInteractionRegistry::
RegisterMe(G4HadronicInteraction * aModel)
{
  if(nModels > 0) {
    for (G4int i=0; i<nModels; i++) {
      if( aModel == allModels[i] ) return;
    }
  }
  //G4cout << "Register model <" << aModel->GetModelName() 
  //<< ">  " << nModels << G4endl;
  allModels.push_back(aModel);
  nModels++;
}

void G4HadronicInteractionRegistry::
RemoveMe(G4HadronicInteraction * aModel)
{
  if(nModels > 0) {
    for (G4int i=0; i<nModels; i++) {
      if( aModel == allModels[i] ) {
	//G4cout << "DeRegister model <" << aModel->GetModelName() 
	//<< ">  " << i << G4endl;
	allModels[i] = 0;
	return;
      }
    }
  }
}
