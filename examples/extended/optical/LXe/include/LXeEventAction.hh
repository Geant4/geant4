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

#ifndef LXeEventAction_h
#define LXeEventAction_h 1

#include "LXeEventMessenger.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Event;
class RecorderBase;

class LXeEventAction : public G4UserEventAction
{
public:
  LXeEventAction(RecorderBase*);
  ~LXeEventAction();
  
public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  
  void SetSaveThreshold(G4int save);

  void SetEventVerbose(G4int v){verbose=v;}

  void SetPMTThreshold(G4int t){pmtThreshold=t;}

  void SetForceDrawPhotons(G4bool b){forcedrawphotons=b;}
  void SetForceDrawNoPhotons(G4bool b){forcenophotons=b;}

private:
  RecorderBase* recorder;
  LXeEventMessenger* eventMessenger;

  G4int              saveThreshold;

  G4int              scintCollID;
  G4int              pmtCollID;

  G4int              verbose;
  
  G4int              pmtThreshold;
  
  G4bool forcedrawphotons;
  G4bool forcenophotons;

};

#endif

    
