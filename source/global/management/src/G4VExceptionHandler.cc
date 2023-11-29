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
// G4VExceptionHandler class implementation
//
// Author: M.Asai, August 2002
// --------------------------------------------------------------------

#include "G4VExceptionHandler.hh"
#include "G4StateManager.hh"

G4VExceptionHandler::G4VExceptionHandler()
{
  G4StateManager* stateManager = G4StateManager::GetStateManager();
  stateManager->SetExceptionHandler(this);
}

G4VExceptionHandler::G4VExceptionHandler(const G4VExceptionHandler& right)
{
  *this = right;
}

G4VExceptionHandler& G4VExceptionHandler::operator=(
  const G4VExceptionHandler& right)
{
  if(&right == this)
  {
    return *this;
  }
  *this = right;
  return *this;
}

G4bool G4VExceptionHandler::operator==(const G4VExceptionHandler& right) const
{
  return (this == &right);
}

G4bool G4VExceptionHandler::operator!=(const G4VExceptionHandler& right) const
{
  return (this != &right);
}
