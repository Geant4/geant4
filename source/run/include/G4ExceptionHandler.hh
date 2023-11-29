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
// G4ExceptionHandler
//
// Class description:
//
// Abstract base class which needs to be notified when G4Exception occurs.
// The concrete class object derived from this class will be automatically
// registered to G4StateManager and the virtual method Notify() will be
// invoked when G4Exception occurs.

// Author: M.Asai - August 2002
// --------------------------------------------------------------------
#ifndef G4ExceptionHandler_hh
#define G4ExceptionHandler_hh 1

#include "G4ExceptionSeverity.hh"
#include "G4VExceptionHandler.hh"
#include "globals.hh"

class G4ExceptionHandler : public G4VExceptionHandler
{
  public:
    G4ExceptionHandler() = default;
    ~G4ExceptionHandler() override = default;
    G4bool operator==(const G4ExceptionHandler& right) const;
    G4bool operator!=(const G4ExceptionHandler& right) const;

    G4ExceptionHandler(const G4ExceptionHandler&) = delete;
    G4ExceptionHandler& operator=(const G4ExceptionHandler&) = delete;

    // Will be invoked by G4StateManager when G4Exception occurs.
    // If TRUE returned, core dump is generated, while if FALSE,
    // the program execution continues.
    G4bool Notify(const char* originOfException, const char* exceptionCode,
                  G4ExceptionSeverity severity, const char* description) override;

  private:
    void DumpTrackInfo();
};

#endif
