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
// $Id: G4UIsession.hh,v 1.4 2001-07-11 10:01:14 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// $id$

#ifndef G4UIsession_h
#define G4UIsession_h 1

#include "G4coutDestination.hh"
#include "globals.hh"

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
      virtual ~G4UIsession();

      virtual G4UIsession * SessionStart();
      // This method will be invoked by main().
      // Optionally, it can be invoked by another session.
      
      virtual void PauseSessionStart(G4String Prompt);
      // This method will be invoked by G4UImanager
      // when G4kernel becomes to Pause state.
      
      virtual G4int ReceiveG4cout(G4String coutString);
      virtual G4int ReceiveG4cerr(G4String cerrString);
      // These two methods will be invoked by G4strstreambuf.

};



#endif

