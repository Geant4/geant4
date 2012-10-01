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
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		SSAEventAction.hh
//
// Version:		0.b.4
// Date:		16/08/99
// Author:		F Lei
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// DESCRIPTION
// -----------
//
// The function of the SSAEventAction is to instantiate the SSAEventMessenger
// class and (at the EndOfEventAction) to draw the event trajectories according 
// to the drawFlag status. 
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
// SSAEventAction ()
//    Constructor:  Instatiates the SSAEventMessenger class.
//
// ~SSAEventAction ()
//    Destructor:  Deletes the SSAEventMessenger class.
//
// void BeginOfEventAction (const G4Run *aRun)
//    Do nothing    
//
// void EndOfRunAction (const G4Run *aRun)
//    Draw the event trajectories.
//
// void SetDrawFlag(G4String val) 
//    to set the drawFlag using val.
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 16 August 1999, F Lei, DERA UK
// Adapted from a verson by Bill Lockman, SLAC, to whom all credits go:
//
// $Id: exrdm01EventAction.hh,v 1.2 2006-12-13 15:46:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#ifndef exrdm01EventAction_h
#define exrdm01EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class exrdm01EventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class exrdm01EventAction : public G4UserEventAction
{
  public:
    exrdm01EventAction();
   ~exrdm01EventAction();

  public:
    void BeginOfEventAction(const G4Event* anEvent);
    void EndOfEventAction(const G4Event* anEvent);
    
    void SetDrawFlag(G4String val)  {drawFlag = val;};
    
  private:
    G4String drawFlag;                         // control the drawing of event
    exrdm01EventActionMessenger*  eventMessenger;
};

#endif

    




