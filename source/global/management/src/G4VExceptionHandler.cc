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
// $Id: G4VExceptionHandler.cc,v 1.1 2002-08-19 18:20:12 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      ---------------- G4VExceptionHandler ----------------
//             by Makoto Asai (August 2002)
// ------------------------------------------------------------

#include "G4VExceptionHandler.hh"
#include "G4StateManager.hh"

G4VExceptionHandler::G4VExceptionHandler() 
{
  G4StateManager * stateManager = G4StateManager::GetStateManager();
  stateManager->SetExceptionHandler(this);
}

G4VExceptionHandler::~G4VExceptionHandler()
{
}

G4VExceptionHandler::G4VExceptionHandler(const G4VExceptionHandler &right)
{
   *this = right;
}

G4VExceptionHandler& G4VExceptionHandler::operator=(const G4VExceptionHandler &right)
{
   if (&right == this) return *this;
   *this = right;
   return *this;
}

G4int G4VExceptionHandler::operator==(const G4VExceptionHandler &right) const
{
   return (this == &right);
}

G4int G4VExceptionHandler::operator!=(const G4VExceptionHandler &right) const
{
   return (this != &right);
}
