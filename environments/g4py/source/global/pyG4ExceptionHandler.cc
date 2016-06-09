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
// $Id: pyG4ExceptionHandler.cc,v 1.1 2006-11-21 05:58:33 kmura Exp $
// $Name: not supported by cvs2svn $
// ====================================================================
//   pyG4ExceptionHandler.cc
//
//                                         2005 Q
// ====================================================================
#include <boost/python.hpp>
#include "G4VExceptionHandler.hh"
#include "G4StateManager.hh"
#include <stdlib.h>

using namespace boost::python;

// ====================================================================
// Python Exception Handler
// ====================================================================
class PyG4ExceptionHandler : public G4VExceptionHandler {
public:
  PyG4ExceptionHandler():G4VExceptionHandler() {
    G4StateManager* sm= G4StateManager::GetStateManager();
    sm-> SetExceptionHandler(this);
  }
  ~PyG4ExceptionHandler() { }

  G4bool Notify(const char* originOfException,
                const char* exceptionCode,
                G4ExceptionSeverity severity,
                const char* description) {

    std::ostringstream message;
    message << "*** G4Exception : " << exceptionCode << G4endl
            << "      issued by : " << originOfException << G4endl
            << description << G4endl;
    
    switch(severity) {
    case FatalException:
      PyErr_SetString(PyExc_AssertionError,
                      "*** Fatal Exception ***");
      PyErr_Print();
      G4cerr << message.str() << G4endl;
      break;

    case FatalErrorInArgument:
      PyErr_SetString(PyExc_ValueError,
                      "*** Fatal Error In Argument ***");
      PyErr_Print();
      G4cerr << message.str() << G4endl;
      break;

    case RunMustBeAborted:
      PyErr_SetString(PyExc_RuntimeError,
                      "*** Run Must Be Aborted ***");
      PyErr_Print();
      G4cerr << message.str() << G4endl;
      break;

    case EventMustBeAborted:
      PyErr_SetString(PyExc_RuntimeError,
                      "*** Event Must Be Aborted ***");
      PyErr_Print();
      G4cerr << message.str() << G4endl;
      break;

    default:
      PyErr_Warn(PyExc_RuntimeWarning,
                 "*** This is just a warning message. ***");
      G4cerr << message.str() << G4endl;
      break;
    }

    // anyway, no abort.
    return false;
  }
};

// ====================================================================
// module definition
// ====================================================================
void export_G4ExceptionHandler()
{
  class_<PyG4ExceptionHandler, boost::noncopyable>
    ("G4ExceptionHandler", "exception handler")
    ;
}
