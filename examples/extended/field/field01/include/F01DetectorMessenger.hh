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
// $Id: F01DetectorMessenger.hh,v 1.3 2001-10-25 08:10:20 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F01DetectorMessenger_h
#define F01DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class F01DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class F01DetectorMessenger: public G4UImessenger
{
  public:
    F01DetectorMessenger(F01DetectorConstruction* );
   ~F01DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    F01DetectorConstruction*   F01Detector;
    
    G4UIdirectory*             F01detDir;

    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsRadCmd;

    G4UIcmdWithADoubleAndUnit* AbsZposCmd;

    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldZCmd;
    G4UIcmdWithADoubleAndUnit* WorldRCmd;

    G4UIcmdWithoutParameter*   UpdateCmd;

};

#endif

