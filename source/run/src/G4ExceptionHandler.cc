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
// $Id: G4ExceptionHandler.cc,v 1.2 2002-12-04 21:52:40 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      ---------------- G4ExceptionHandler ----------------
//             by Makoto Asai (August 2002)
// ------------------------------------------------------------

#include "G4ExceptionHandler.hh"
#include "G4StateManager.hh"
#include "G4RunManager.hh"
#include "G4ios.hh"
#include <stdlib.h>
#include "G4String.hh"

G4ExceptionHandler::G4ExceptionHandler() 
{
}

G4ExceptionHandler::~G4ExceptionHandler()
{
}

G4ExceptionHandler::G4ExceptionHandler(const G4ExceptionHandler &right)
{
   *this = right;
}

G4ExceptionHandler& G4ExceptionHandler::operator=(const G4ExceptionHandler &right)
{
   if (&right == this) return *this;
   *this = right;
   return *this;
}

G4int G4ExceptionHandler::operator==(const G4ExceptionHandler &right) const
{
   return (this == &right);
}

G4int G4ExceptionHandler::operator!=(const G4ExceptionHandler &right) const
{
   return (this != &right);
}

G4bool G4ExceptionHandler::Notify(const char* originOfException,
                        const char* exceptionCode,
                        G4ExceptionSeverity severity,
                        const char* description)
{
  G4cerr << G4endl;
  G4cerr << "*** G4Exception : " << exceptionCode << G4endl;
  G4cerr << "      issued by : " << originOfException << G4endl;
  G4cerr << description << G4endl;
  G4bool abortionForCoreDump = false;
  G4ApplicationState aps = G4StateManager::GetStateManager()->GetCurrentState();
  switch(severity)
  {
   case FatalException:
    G4cerr << "*** Fatal Exception *** core dump ***";
    abortionForCoreDump = true;
    break;
   case FatalErrorInArgument:
    G4cerr << "*** Fatal Error In Argument *** core dump ***";
    abortionForCoreDump = true;
    break;
   case RunMustBeAborted:
    if(aps==G4State_GeomClosed || aps==G4State_EventProc)
    {
      G4cerr << "*** Run Must Be Aborted ";
      G4RunManager::GetRunManager()->AbortRun(false);
    }
    abortionForCoreDump = false;
    break;
   case EventMustBeAborted:
    if(aps==G4State_EventProc)
    {
      G4cerr << "*** Event Must Be Aborted ";
      G4RunManager::GetRunManager()->AbortEvent();
    }
    abortionForCoreDump = false;
    break;
   default:
    G4cerr << "*** This is just a warning message.";
    abortionForCoreDump = false;
    break;
  }
  G4cerr << G4endl;
  return abortionForCoreDump;
}


