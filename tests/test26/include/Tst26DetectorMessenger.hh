//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Tst26DetectorMessenger.hh,v 1.4 2006-06-29 21:53:07 gunter Exp $
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

