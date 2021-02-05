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

#include "G4InteractorMessenger.hh"

#include "G4VInteractiveSession.hh"

/***************************************************************************/
G4VInteractiveSession::G4VInteractiveSession ()
{
  messenger = new G4InteractorMessenger(this);
}

/***************************************************************************/
G4VInteractiveSession::~G4VInteractiveSession() 
{
  delete messenger;
}

/***************************************************************************/
void G4VInteractiveSession::AddMenu (const char*,const char*)
{
}

/***************************************************************************/
void G4VInteractiveSession::AddButton (const char*,const char*,const char*)
{
}

/***************************************************************************/
void G4VInteractiveSession::DefaultIcons (bool)
{
}

/***************************************************************************/
void G4VInteractiveSession::AddIcon (const char*,const char*,const char*,const char*)
{
}

/***************************************************************************/
void G4VInteractiveSession::OutputStyle(const char*,const char*,const char*)
{
}

/***************************************************************************/
void G4VInteractiveSession::AddInteractor (G4String a_name,
                                           G4Interactor a_interactor)
{
  interactors[a_name] = a_interactor;
}

/***************************************************************************/
G4Interactor G4VInteractiveSession::GetInteractor (G4String a_name)
{
  G4interactor_map::iterator it;
  if((it=interactors.find(a_name))==interactors.end()) return NULL;
  return (*it).second;  
}
