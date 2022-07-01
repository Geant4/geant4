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
// G4UIsession
//
// Class description:
//
// This is a base class of all (G)UI sessions.
// SessionStart() method should be called to start the session

// Author: Makoto Asai, 1998
// --------------------------------------------------------------------
#ifndef G4UIsession_hh
#define G4UIsession_hh 1

#include "G4coutDestination.hh"
#include "globals.hh"
#include "icomsdefs.hh"

class G4UIsession : public G4coutDestination
{
  public:

    G4UIsession();
    G4UIsession(G4int iBatch);
    ~G4UIsession() override;

    virtual G4UIsession* SessionStart();
      // This method will be invoked by main().
      // Optionally, it can be invoked by another session

    virtual void PauseSessionStart(const G4String& Prompt);
      // This method will be invoked by G4UImanager
      // when the kernel comes to Pause state

    G4int ReceiveG4cout(const G4String& coutString) override;
    G4int ReceiveG4cerr(const G4String& cerrString) override;
    // These two methods will be invoked by G4strstreambuf

    static G4int InSession();
    inline G4int GetLastReturnCode() const { return lastRC; }

  protected:

    G4ICOMS_DLL static G4int inSession;
    G4int ifBatch = 0;
    G4int lastRC = 0;
};

#endif
