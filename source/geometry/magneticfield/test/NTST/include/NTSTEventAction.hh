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


// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: NTSTEventAction.hh,v 1.2 2003-12-09 15:35:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTEventAction_h
#define NTSTEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Timer;
class G4Event;
class NTSTEventActionMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTEventAction : public G4UserEventAction
{
public:
  NTSTEventAction();
  ~NTSTEventAction();
  
public:
  void BeginOfEventAction(const G4Event* aEvent);
  void EndOfEventAction(const G4Event* aEvent);
  
  void SetDrawFlag(G4String val)  {drawFlag = val;};
  
private:
  G4Timer* EventTime;
  double MeanUserEventTime;
  double RmsUserEventTime;
  double MeanRealEventTime;
  double RmsRealEventTime;
  double NumberOfEvents;
  double MeanVertices;
  double RmsVertices;
  double MeanTracks;
  double RmsTracks;
  G4int    calorimeterCollID;                // Hits collection ID
  
  G4String drawFlag;                         // control the drawing of event
  NTSTEventActionMessenger*  eventMessenger;
};

#endif

    

