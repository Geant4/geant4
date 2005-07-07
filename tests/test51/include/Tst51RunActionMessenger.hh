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
// $Id: Tst51RunActionMessenger.hh
// GEANT4 tag $Name: xray_fluo-V04-01-03
//
//
// History:
// -----------
//
// -------------------------------------------------------------------


#ifndef Tst51RunActionMessenger_h
#define Tst51RunActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst51RunAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class Tst51RunActionMessenger: public G4UImessenger
{
  
public:

  Tst51RunActionMessenger(Tst51RunAction*);
  ~Tst51RunActionMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:

  Tst51RunAction*          fTst51RunAction;

  G4UIdirectory* analysisDir;
  G4UIcmdWithAString*  fileNameCmd;
};

#endif








