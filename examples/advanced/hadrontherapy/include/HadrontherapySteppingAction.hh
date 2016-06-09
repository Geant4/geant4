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
// $Id: HadrontherapySteppingAction.cc,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
#ifndef HadrontherapySteppingAction_h
#define HadrontherapySteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "globals.hh"
#include <iostream.h>

class HadrontherapyDetectorConstruction;
class HadrontherapyRunAction;
class HadrontherapyEventAction;
class HadrontherapySteppingMessenger;

// ---------------------------------------------------------------
#ifndef G4NOHIST
 class HepTupleManager;
 class HepHistogram;
#endif

class HadrontherapySteppingAction : public G4UserSteppingAction
{
public:
  HadrontherapySteppingAction( HadrontherapyEventAction* );
  ~HadrontherapySteppingAction();
  
  void UserSteppingAction(const G4Step*);

  G4int event_id;
  G4int Controllo;

private:
  HadrontherapyDetectorConstruction* detector;
  HadrontherapyEventAction*          eventaction;
  HadrontherapyRunAction*            runaction;
  HadrontherapySteppingMessenger*    steppingMessenger;   
};
#endif
