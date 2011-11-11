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
// $Id: pyG4Exception.cc,v 1.1 2006-11-20 05:57:16 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4Exception.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4StateManager.hh"
#include "G4ExceptionSeverity.hh"

using namespace boost::python;

// ====================================================================
// thin wrappers
// ====================================================================
namespace pyG4Exception {

//////////////////////////////////////////////////
void f2_G4Exception(const char* originOfException,
                    const char* exceptionCode,
                    G4ExceptionSeverity severity,
                    const char* description)
//////////////////////////////////////////////////
{
  G4VExceptionHandler* exceptionHandler
    = G4StateManager::GetStateManager()-> GetExceptionHandler();
  G4bool toBeAborted = true;
  if(exceptionHandler) {
    toBeAborted = exceptionHandler
      -> Notify(originOfException,exceptionCode,severity,description);
  } else {
    G4String e_banner = "\n!!!!! - !!!!! - !!!!! - !!!!! - !!!!! - !!!!!\n";
    G4String w_banner = "\nwwwww - wwwww - wwwww - wwwww - wwwww - wwwww\n";
    std::ostringstream message;
    message << "\n*** ExceptionHandler is not defined ***\n"
            << "*** G4Exception : " << exceptionCode << G4endl
            << "      issued by : " << originOfException << G4endl
            << description << G4endl;
    switch(severity) {
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
      G4cout << w_banner << message.str() 
             << "*** This is just a warning message. ***"
             << w_banner;
      toBeAborted = false;
      break;
    }
  }

  if(toBeAborted) {
    if(G4StateManager::GetStateManager()->SetNewState(G4State_Abort)) {
      G4cerr << G4endl << "*** G4Exception: Aborting execution ***" << G4endl;
      PyErr_SetString(PyExc_RuntimeError, description);
      PyErr_Print();
    } else {
      G4cerr << G4endl << "*** G4Exception: Abortion suppressed ***"
             << G4endl << "*** No guarantee for further execution ***" 
             << G4endl;
    }
  }
}

}

using namespace pyG4Exception;

// ====================================================================
// module definition
// ====================================================================
void export_G4Exception()
{
  def("G4Exception",   f2_G4Exception);
}
