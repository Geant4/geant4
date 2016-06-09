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
// $Id: MedLinacDetectorMessenger.hh,v 1.3 2005/07/03 23:27:36 mpiergen Exp $
//
//
// Code developed by: M. Piergentili

#ifndef MedLinacDetectorMessenger_h
#define MedLinacDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class MedLinacDetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//****************************************************************************

class MedLinacDetectorMessenger: public G4UImessenger
{
  public:
    MedLinacDetectorMessenger(MedLinacDetectorConstruction*);
   ~MedLinacDetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    MedLinacDetectorConstruction* MedLinacDetector;

    G4UIdirectory*             MedLinacDir;
    G4UIdirectory*             X1Dir;
    G4UIdirectory*             X2Dir;
    G4UIdirectory*             Y1Dir;
    G4UIdirectory*             Y2Dir;
    G4UIcmdWithADoubleAndUnit* JawX1PosCmd;
    G4UIcmdWithADoubleAndUnit* JawX2PosCmd;
    G4UIcmdWithADoubleAndUnit* JawY1PosCmd;
    G4UIcmdWithADoubleAndUnit* JawY2PosCmd;

    G4UIcmdWithADoubleAndUnit* PhantomDimCmd;
    G4UIcmdWithAnInteger* NVoxelsCmd;
    G4UIcmdWithADoubleAndUnit* MaxStepCmd;
   
    G4UIcmdWithoutParameter*   UpdateCmd;
};

//****************************************************************************
#endif

