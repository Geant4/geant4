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
// $Id: HadrontherapyEventAction.hh,v 1.0
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------

#ifndef HadrontherapyEventAction_h
#define HadrontherapyEventAction_h 1
#include "G4UserEventAction.hh"
#include "globals.hh"

class HadrontherapyRunAction;
class HadrontherapyEventActionMessenger;
class HadrontherapyHit;
// ---------------------------------------------------------------
class HadrontherapyEventAction : public G4UserEventAction
{
public:
  HadrontherapyEventAction(HadrontherapyRunAction* runAction );
  ~HadrontherapyEventAction();

public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  void setEventVerbose(G4int level);
  void CountStepsCharged() ;
  void CountStepsNeutral() ;
  void AddCharged() ;
  void AddNeutral() ;
  void AddE();
  void AddP();   
  void SetTr();
  void SetRef();
  
  G4int GetEventno();
  G4int Trasporto();

private:
  G4int event_id;
  G4int    calorimeterCollID;
  HadrontherapyEventActionMessenger*  eventMessenger;
  HadrontherapyRunAction* p_Run;
  G4int verboselevel;
  G4double nstep;
  G4double energyDep[50000];
};
#endif
