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
/// \file ParRunManager.hh
/// \brief Definition of the ParRunManager class
//
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
