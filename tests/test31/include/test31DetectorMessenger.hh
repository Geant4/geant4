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
#ifndef test31DetectorMessenger_h
#define test31DetectorMessenger_h 1

// -------------------------------------------------------------
//
//
// -------------------------------------------------------------
//      GEANT4 test31
//
//      History: based on object model of
//      2nd December 1995, G.Cosmo
//      ---------- test31DetectorMessenger -------
//              
//  Modified: 05.04.01 Vladimir Ivanchenko new design of test31 
// 
// -------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "globals.hh"
#include "G4UImessenger.hh"

class test31DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class test31DetectorMessenger: public G4UImessenger
{
public: // Without description

    test31DetectorMessenger(test31DetectorConstruction* );
   ~test31DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    test31DetectorConstruction* hDet;
    
    G4UIdirectory*             test31detDir;
    G4UIdirectory*             test31detDir1;
    G4UIdirectory*             test31detDir2;

    G4UIcmdWithAnInteger*      NumOfAbsCmd;
    G4UIcmdWithAString*        AbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* AbsThickCmd;
    G4UIcmdWithADoubleAndUnit* AbsGapCmd;
    G4UIcmdWithADoubleAndUnit* AbsSizYZCmd;
    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldXCmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
    G4UIcmdWithADoubleAndUnit* XMagFieldCmd;
    G4UIcmdWithADoubleAndUnit* YMagFieldCmd;
    G4UIcmdWithADoubleAndUnit* ZMagFieldCmd;
    G4UIcmdWithAString*        HistoCmd;
    G4UIcmdWithAnInteger*      NumOfEvt;
    G4UIcmdWithAnInteger*      verbCmd;
    G4UIcmdWithAnInteger*      intCmd;
    G4UIcmdWithAnInteger*      nhistCmd;
    G4UIcmdWithAnInteger*      nDebugSCmd;
    G4UIcmdWithAnInteger*      nDebugECmd;
    G4UIcmdWithADoubleAndUnit* DeltaECmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif





