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

#ifndef LXeEventMessenger_h
#define LXeEventMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class LXeEventAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithABool;

class LXeEventMessenger: public G4UImessenger
{
public:
  LXeEventMessenger(LXeEventAction*);
  ~LXeEventMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  LXeEventAction*        LXeEvent;
  G4UIcmdWithAnInteger*  saveThresholdCmd;
  G4UIcmdWithAnInteger*  verboseCmd;
  G4UIcmdWithAnInteger*  pmtThresholdCmd;
  G4UIcmdWithABool*      forceDrawPhotonsCmd;
  G4UIcmdWithABool*      forceDrawNoPhotonsCmd;
};

#endif

