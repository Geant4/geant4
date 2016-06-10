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
// $Id: G4VDigitizerModule.cc 80987 2014-05-19 10:50:22Z gcosmo $
//

#include "G4VDigitizerModule.hh"
#include "G4VDigiCollection.hh"
#include "G4DigiManager.hh"

G4VDigitizerModule::G4VDigitizerModule(G4String modName)
:verboseLevel(0)
{
  moduleName = modName;
  DigiManager = G4DigiManager::GetDMpointer();
}

G4VDigitizerModule::~G4VDigitizerModule()
{;}

int G4VDigitizerModule::operator==(const G4VDigitizerModule &right) const
{ return (moduleName==right.moduleName); }

int G4VDigitizerModule::operator!=(const G4VDigitizerModule &right) const
{ return (moduleName!=right.moduleName); }

void G4VDigitizerModule::StoreDigiCollection(G4VDigiCollection* aDC)
{
  G4String DCnam = moduleName;
  DCnam += "/";
  DCnam += aDC->GetName();
  G4int DCID = DigiManager->GetDigiCollectionID(DCnam);
  if(DCID>=0) StoreDigiCollection(DCID,aDC);
}

void G4VDigitizerModule::StoreDigiCollection(G4int DCID,G4VDigiCollection* aDC)
{
  DigiManager->SetDigiCollection(DCID,aDC);
}


