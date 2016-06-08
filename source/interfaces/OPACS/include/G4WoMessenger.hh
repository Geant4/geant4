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
// $Id: G4WoMessenger.hh,v 1.3.4.1 2001/06/28 19:10:20 gunter Exp $
// GEANT4 tag $Name:  $
//
#ifndef G4WoMessenger_h
#define G4WoMessenger_h 1

#ifdef G4UI_BUILD_WO_SESSION

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

#include <OShell.h>

class G4WoMessenger: public G4UImessenger
{
  public:
    G4WoMessenger(OShell);
    ~G4WoMessenger();
  public:
    void     SetNewValue    (G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);
  private:
    OShell shell;
  private: //commands
    G4UIdirectory*           WoDirectory;
    G4UIcmdWithoutParameter* sendExitCmd;
    G4UIcmdWithAString*      oshCmd;
};

#endif

#endif





