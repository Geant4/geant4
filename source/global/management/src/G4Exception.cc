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
// $Id: G4Exception.cc,v 1.21 2007-11-13 17:35:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
// G4Exception
//
// Global error function prints string to G4cerr (or G4cout in case of
// warning). May abort program according to severity.
// ----------------------------------------------------------------------

#include "G4ios.hh"
#include "G4String.hh"
#include "G4StateManager.hh"

void G4Exception(const char* emessage)
{
   if(emessage)  { G4cerr << emessage << G4endl; }
   if(G4StateManager::GetStateManager()->SetNewState(G4State_Abort,emessage))
   {
     G4cerr << G4endl << "*** G4Exception: Aborting execution ***" << G4endl;
     abort();
   }
   else
   {
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
    G4String e_banner = "\n!!!!! - !!!!! - !!!!! - !!!!! - !!!!! - !!!!!\n";
    G4String w_banner = "\nwwwww - wwwww - wwwww - wwwww - wwwww - wwwww\n";
    std::ostringstream message;
    message << "\n*** ExceptionHandler is not defined ***\n"
            << "*** G4Exception : " << exceptionCode << G4endl
            << "      issued by : " << originOfException << G4endl
            << description << G4endl;
    switch(severity)
    {
     case FatalException:
      G4cerr << e_banner << message.str() << "*** Fatal Exception ***"
             << e_banner;
      break;
     case FatalErrorInArgument:
      G4cerr << e_banner << message.str() << "*** Fatal Error In Argument ***"
             << e_banner;
      break;
     case RunMustBeAborted:
      G4cerr << e_banner << message.str() << "*** Run Must Be Aborted ***"
             << e_banner;
      break;
     case EventMustBeAborted:
      G4cerr << e_banner << message.str() << "*** Event Must Be Aborted ***"
             << e_banner;
      break;
     default:
      G4cout << w_banner << message.str() << "*** This is just a warning message. ***"
             << w_banner;
      toBeAborted = false;
      break;
    }
  }
  if(toBeAborted)
  {
   if(G4StateManager::GetStateManager()->SetNewState(G4State_Abort))
   {
     G4cerr << G4endl << "*** G4Exception: Aborting execution ***" << G4endl;
     abort();
   }
   else
   {
     G4cerr << G4endl << "*** G4Exception: Abortion suppressed ***"
            << G4endl << "*** No guarantee for further execution ***" << G4endl;
   }
  }
}

void G4Exception(std::string emessage)
{
  G4Exception(emessage.c_str());
}

void G4Exception(G4String emessage)
{
  G4Exception(emessage.c_str());
}

//typedef std::ostringstream G4ExceptionDescription;

void G4Exception(const char* originOfException,
                 const char* exceptionCode,
                 G4ExceptionSeverity severity,
                 G4ExceptionDescription & description)
{
  G4String des = description.str();
  G4Exception(originOfException, exceptionCode, severity, des.c_str());
}

void G4Exception(const char* originOfException,
                 const char* exceptionCode,
                 G4ExceptionSeverity severity,
                 G4ExceptionDescription & description,
                 const char* comments)
{
  description << comments << G4endl;
  G4Exception(originOfException, exceptionCode, severity, description);
}
