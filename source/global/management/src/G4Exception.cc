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
// G4Exception implementation
//
// Authors: G.Cosmo, M.Asai - May 1999 - First implementation
// --------------------------------------------------------------------

#include "G4Exception.hh"
#include "G4ios.hh"
#include "G4StateManager.hh"
#include "G4VExceptionHandler.hh"

namespace {

inline const G4String G4ExceptionErrBannerStart()
{
  return "\n-------- EEEE ------- G4Exception-START -------- EEEE -------\n";
}
inline const G4String G4ExceptionWarnBannerStart()
{
  return "\n-------- WWWW ------- G4Exception-START -------- WWWW -------\n";
}

inline const G4String G4ExceptionErrBannerEnd()
{
  return "\n-------- EEEE ------- G4Exception-END -------- EEEE -------\n";
}
inline const G4String G4ExceptionWarnBannerEnd()
{
  return "\n-------- WWWW ------- G4Exception-END -------- WWWW -------\n";
}

} // namespace

// --------------------------------------------------------------------
void G4Exception(const char* originOfException, const char* exceptionCode,
                 G4ExceptionSeverity severity, const char* description)
{
  G4VExceptionHandler* exceptionHandler =
    G4StateManager::GetStateManager()->GetExceptionHandler();
  G4bool toBeAborted = true;
  if(exceptionHandler != nullptr)
  {
    toBeAborted = exceptionHandler->Notify(originOfException, exceptionCode,
                                           severity, description);
  }
  else
  {
    static const G4String& es_banner = G4ExceptionErrBannerStart();
    static const G4String& ee_banner = G4ExceptionErrBannerEnd();
    static const G4String& ws_banner = G4ExceptionWarnBannerStart();
    static const G4String& we_banner = G4ExceptionWarnBannerEnd();
    std::ostringstream message;
    message << "\n*** ExceptionHandler is not defined ***\n"
            << "*** G4Exception : " << exceptionCode << G4endl
            << "      issued by : " << originOfException << G4endl
            << description << G4endl;
    switch(severity)
    {
      case FatalException:
        G4cerr << es_banner << message.str() << "*** Fatal Exception ***"
               << ee_banner << G4endl;
        break;
      case FatalErrorInArgument:
        G4cerr << es_banner << message.str()
               << "*** Fatal Error In Argument ***" << ee_banner << G4endl;
        break;
      case RunMustBeAborted:
        G4cerr << es_banner << message.str() << "*** Run Must Be Aborted ***"
               << ee_banner << G4endl;
        break;
      case EventMustBeAborted:
        G4cerr << es_banner << message.str() << "*** Event Must Be Aborted ***"
               << ee_banner << G4endl;
        break;
      default:
        G4cout << ws_banner << message.str()
               << "*** This is just a warning message. ***" << we_banner
               << G4endl;
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
      G4cerr << G4endl << "*** G4Exception: Abortion suppressed ***" << G4endl
             << "*** No guarantee for further execution ***" << G4endl;
    }
  }
}

// --------------------------------------------------------------------
void G4Exception(const char* originOfException, const char* exceptionCode,
                 G4ExceptionSeverity severity,
                 G4ExceptionDescription& description)
{
  G4String des = description.str();
  G4Exception(originOfException, exceptionCode, severity, des.c_str());
}

// --------------------------------------------------------------------
void G4Exception(const char* originOfException, const char* exceptionCode,
                 G4ExceptionSeverity severity,
                 G4ExceptionDescription& description, const char* comments)
{
  description << comments << G4endl;
  G4Exception(originOfException, exceptionCode, severity, description);
}
