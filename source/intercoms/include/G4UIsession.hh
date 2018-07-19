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
// $Id: G4UIsession.hh 106172 2017-09-15 13:03:57Z gcosmo $
//
// $id$

#ifndef G4UIsession_h
#define G4UIsession_h 1

#include "G4coutDestination.hh"
#include "globals.hh"
#include "icomsdefs.hh"

// class description:
//
//  This is a base class of all (G)UI session.
//  SessionStart() method should be called to start the session.
//

class G4UIsession : public G4coutDestination
{
  // Base class of UI/GUI session
  
  public:
      G4UIsession();
      G4UIsession(G4int iBatch);
      virtual ~G4UIsession();

      virtual G4UIsession * SessionStart();
      // This method will be invoked by main().
      // Optionally, it can be invoked by another session.
      
      virtual void PauseSessionStart(const G4String& Prompt);
      // This method will be invoked by G4UImanager
      // when G4kernel becomes to Pause state.
      
      virtual G4int ReceiveG4cout(const G4String& coutString);
      virtual G4int ReceiveG4cerr(const G4String& cerrString);
      // These two methods will be invoked by G4strstreambuf.

  protected:
      G4ICOMS_DLL static G4int inSession;
      G4int ifBatch;
  public:
      static G4int InSession() { return inSession; }
};



#endif

