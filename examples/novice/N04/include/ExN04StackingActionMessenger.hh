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

#ifndef ExN04StackingActionMessenger_h
#define ExN04StackingActionMessenger_h 1

class ExN04StackingAction;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;

#include "G4UImessenger.hh"
#include "globals.hh"

class ExN04StackingActionMessenger: public G4UImessenger
{
  public:
    ExN04StackingActionMessenger(ExN04StackingAction* msa);
    ~ExN04StackingActionMessenger();
    
  public:
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    ExN04StackingAction * myAction;
    
  private: //commands
    G4UIcmdWithAnInteger * muonCmd;
    G4UIcmdWithAnInteger * isomuonCmd;
    G4UIcmdWithAnInteger * isoCmd;
    G4UIcmdWithADoubleAndUnit * roiCmd;
    
};

#endif


