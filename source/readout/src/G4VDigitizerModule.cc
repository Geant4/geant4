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
// $Id: G4VDigitizerModule.cc,v 1.3 2001-07-11 10:08:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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


