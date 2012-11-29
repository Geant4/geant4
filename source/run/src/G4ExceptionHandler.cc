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
// $Id$
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

G4ExceptionHandler::G4ExceptionHandler(const G4ExceptionHandler &)
:G4VExceptionHandler()
{
}

G4ExceptionHandler& G4ExceptionHandler::operator=(const G4ExceptionHandler &)
{
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
  static const G4String es_banner
    = "\n-------- EEEE ------- G4Exception-START -------- EEEE -------\n";
  static const G4String ee_banner
    = "\n-------- EEEE -------- G4Exception-END --------- EEEE -------\n";
  static const G4String ws_banner
    = "\n-------- WWWW ------- G4Exception-START -------- WWWW -------\n";
  static const G4String we_banner
    = "\n-------- WWWW -------- G4Exception-END --------- WWWW -------\n";
  std::ostringstream message;
  message << "*** G4Exception : " << exceptionCode << G4endl
          << "      issued by : " << originOfException << G4endl
          << description << G4endl;
  G4bool abortionForCoreDump = false;
  G4ApplicationState aps = G4StateManager::GetStateManager()->GetCurrentState();
  switch(severity)
  {
   case FatalException:
    G4cerr << es_banner << message.str() << "*** Fatal Exception *** core dump ***"
           << ee_banner << G4endl;
    abortionForCoreDump = true;
    break;
   case FatalErrorInArgument:
    G4cerr << es_banner << message.str() << "*** Fatal Error In Argument *** core dump ***"
           << ee_banner << G4endl;
    abortionForCoreDump = true;
    break;
   case RunMustBeAborted:
    if(aps==G4State_GeomClosed || aps==G4State_EventProc)
    {
      G4cerr << es_banner << message.str() << "*** Run Must Be Aborted ***"
             << ee_banner << G4endl;
      G4RunManager::GetRunManager()->AbortRun(false);
    }
    abortionForCoreDump = false;
    break;
   case EventMustBeAborted:
    if(aps==G4State_EventProc)
    {
      G4cerr << es_banner << message.str() << "*** Event Must Be Aborted ***"
             << ee_banner << G4endl;
      G4RunManager::GetRunManager()->AbortEvent();
    }
    abortionForCoreDump = false;
    break;
   default:
    G4cout << ws_banner << message.str() << "*** This is just a warning message. ***"
           << we_banner << G4endl;
    abortionForCoreDump = false;
    break;
  }
  return abortionForCoreDump;
}
