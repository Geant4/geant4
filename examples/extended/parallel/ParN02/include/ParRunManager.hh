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
// $Id: ParRunManager.hh,v 1.2 2002-06-06 17:05:02 cooperma Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------

#ifndef ParRunManager_h
#define ParRunManager_h 1

#include "G4RunManager.hh"
#include "G4StateManager.hh"
#include "G4UserEventAction.hh"
#include "topc.h"

class ParRunManager : public G4RunManager
{

  protected:
    virtual void DoEventLoop(G4int n_event,
                             const char* macroFile=0,
                             G4int n_select=-1);

  private:  // TOP-C callbacks
    TOPC_BUF GenerateEventInput();
    TOPC_BUF DoEvent( void * input_buf );
    TOPC_ACTION CheckEventResult( void * input_buf, void * ignore );

    static ParRunManager* myRunManager;
    static TOPC_BUF MyGenerateEventInput();
    static TOPC_BUF MyDoEvent( void * input_buf );
    static TOPC_ACTION MyCheckEventResult( void * input_buf, void * ignore );

  // Copies of local variables of G4RunManager::DoEventLoop()
    static G4StateManager* stateManager;
    static G4int n_event;
    static G4int n_select;
    static G4String msg;

    G4int i_event;  // Used in lieu of i_event in RunManager::DoEventLoop()

    inline void ImportDoEventLoopLocals(
	 G4StateManager* stateManager, G4int n_event,
	 G4int n_select, G4String msg )
      {
        ParRunManager::stateManager = stateManager;
        ParRunManager::n_event = n_event;
        ParRunManager::n_select = n_select;
        ParRunManager::msg = msg;
        ParRunManager::myRunManager = this;
      }

  // Copy of userEventAction:  see ParRunManager.cc for further comments
    G4UserEventAction * origUserEventAction;

};

#endif
