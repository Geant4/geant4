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
// $Id: NTSTPhysicsListMessenger.hh,v 1.2 2003-12-09 15:35:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef NTSTPhysicsListMessenger_h
#define NTSTPhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class NTSTPhysicsList;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithADoubleAndUnit;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class NTSTPhysicsListMessenger: public G4UImessenger
{
public:
  NTSTPhysicsListMessenger(NTSTPhysicsList*);
  ~NTSTPhysicsListMessenger();
  
  void SetNewValue(G4UIcommand*, G4String);
  
private:
  NTSTPhysicsList*           NTSTList;
  G4UIcmdWithoutParameter*   ProcessCmd;
  G4UIcmdWithABool*          useBgsTranCmd;
  G4UIcmdWithADoubleAndUnit* MinimumEnergyCutCmd;
  G4UIcmdWithADoubleAndUnit* MaximumEnergyCutCmd;
  G4UIcmdWithADoubleAndUnit* CutCmd;
  G4UIcmdWithADoubleAndUnit* LooperCutCmd;
  
};

#endif

