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
// $Id: Tst26DetectorMessenger.hh,v 1.3 2003-02-06 11:53:27 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst26DetectorMessenger_h
#define Tst26DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Tst26DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWith3Vector;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst26DetectorMessenger: public G4UImessenger
{
  public:
    Tst26DetectorMessenger(Tst26DetectorConstruction* );
   ~Tst26DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:
    Tst26DetectorConstruction* Tst26Detector;
    
    G4UIdirectory*             testemDir;
    G4UIcmdWithAString*        MaterCmd;
    G4UIcmdWithAString*        LBinCmd;
    G4UIcmdWithADoubleAndUnit* l1Cmd;
    G4UIcmdWithADoubleAndUnit* l2Cmd;
    G4UIcmdWithADoubleAndUnit* l3Cmd;
    G4UIcmdWithADoubleAndUnit* l4Cmd;
    G4UIcmdWithADoubleAndUnit* l5Cmd;
    G4UIcmdWithADoubleAndUnit* l6Cmd;
    G4UIcmdWithoutParameter*   UpdateCmd;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

