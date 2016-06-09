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
// $Id: G4Exception.cc,v 1.17 2003/10/31 17:44:47 asaim Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// 
// ----------------------------------------------------------------------
// G4Exception
//
// Global error function prints string to G4cerr, and aborts
// program
//
// History:
// 30.06.95 P.Kent

#include "G4ios.hh"
#include <stdlib.h>
#include "G4String.hh"
#include "G4StateManager.hh"

void G4Exception(const char* s)
{
   if(s)
	{
	    G4cerr << s << G4endl;
	}
   if(G4StateManager::GetStateManager()->SetNewState(G4State_Abort,s)) {
     G4cerr << G4endl << "*** G4Exception: Aborting execution ***" << G4endl;
     abort();
   } else {
     G4cerr << G4endl << "*** G4Exception: Abortion suppressed ***"
            << G4endl << "*** No guarantee for further execution ***" << G4endl;
   }
}

void G4Exception(const char* originOfException,
                        const char* exceptionCode,
                        G4ExceptionSeverity severity,
                        const char* description)
{
  G4VExceptionHandler* exceptionHandler
    = G4StateManager::GetStateManager()->GetExceptionHandler();
  G4bool toBeAborted = true;
  if(exceptionHandler)
  {
    toBeAborted = exceptionHandler
     ->Notify(originOfException,exceptionCode,severity,description);
  }
  else
  {
    G4cerr << G4endl
       << "*** ExceptionHandler is not defined ***" << G4endl;
    G4cerr << G4endl;
    G4cerr << "*** G4Exception : " << exceptionCode << G4endl;
    G4cerr << "      issued by : " << originOfException << G4endl;
    G4cerr << description << G4endl;
    G4cerr << G4endl << "Severity : ";
    switch(severity)
    {
     case FatalException:
      G4cerr << "*** Fatal Exception ***";
      break;
     case FatalErrorInArgument:
      G4cerr << "*** Fatal Error In Argument ***";
      break;
     case RunMustBeAborted:
      G4cerr << "*** Run Must Be Aborted ***";
      break;
     case EventMustBeAborted:
      G4cerr << "*** Event Must Be Aborted ***";
      break;
     default:
      G4cerr << "*** This is just a warning message. ***";
      toBeAborted = false;
      break;
    }
    G4cerr << G4endl;
  }
  if(toBeAborted)
  {
   if(G4StateManager::GetStateManager()->SetNewState(G4State_Abort)) {
     G4cerr << G4endl << "*** G4Exception: Aborting execution ***" << G4endl;
     abort();
   } else {
     G4cerr << G4endl << "*** G4Exception: Abortion suppressed ***"
            << G4endl << "*** No guarantee for further execution ***" << G4endl;
   }
  }
}

void G4Exception(std::string s)
{
  G4Exception(s.c_str());
}

void G4Exception(G4String s)
{
  G4Exception(s.c_str());
}
