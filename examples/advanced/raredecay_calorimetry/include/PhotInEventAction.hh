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
// $Id: PhotInEventAction.hh,v 1.3 2006/06/29 16:24:45 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//

#ifndef PhotInEventAction_h
#define PhotInEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "PhotInStackingAction.hh"
#include "PhotInConstants.hh"
#include "PhotInCalorHit.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

class PhotInEventAction : public G4UserEventAction
{
public:
  PhotInEventAction();
  virtual ~PhotInEventAction();

  virtual void   BeginOfEventAction(const G4Event*);
  virtual void   EndOfEventAction(const G4Event*);
    
  static void SetVerboseLevel(G4int i) { verboseLevel = i; }
  static G4int GetVerboseLevel()       { return verboseLevel; }

private:
  G4int          calorimeterCollID[PhotInDiNSections]; //Collections
  static G4int   verboseLevel;
};


#endif

    
