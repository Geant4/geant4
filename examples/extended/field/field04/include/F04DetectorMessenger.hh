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
// $Id: F04DetectorMessenger.hh,v 1.1 2007-10-30 02:01:47 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef F04DetectorMessenger_h
#define F04DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "F04DetectorConstruction.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class F04DetectorMessenger : public G4UImessenger
{
  public:

    F04DetectorMessenger(F04DetectorConstruction* );
    ~F04DetectorMessenger();
    
    void SetNewValue(G4UIcommand*, G4String);
    
  private:

    F04DetectorConstruction*   Detector;
    
    G4UIdirectory*          detDir;

    G4UIcmdWithAString*        WorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* WorldRCmd;
    G4UIcmdWithADoubleAndUnit* WorldZCmd;

    G4UIcmdWithADoubleAndUnit* CaptureRCmd;
    G4UIcmdWithADoubleAndUnit* CaptureZCmd;
    G4UIcmdWithADoubleAndUnit* CaptureB1Cmd;
    G4UIcmdWithADoubleAndUnit* CaptureB2Cmd;

    G4UIcmdWithADoubleAndUnit* TransferRCmd;
    G4UIcmdWithADoubleAndUnit* TransferZCmd;
    G4UIcmdWithADoubleAndUnit* TransferBCmd;
    G4UIcmdWithADoubleAndUnit* TransferPCmd;

    G4UIcmdWithAString*        TgtMaterCmd;
    G4UIcmdWithADoubleAndUnit* TgtRadCmd;
    G4UIcmdWithADoubleAndUnit* TgtThickCmd;
    G4UIcmdWithADoubleAndUnit* TgtPosCmd;
    G4UIcmdWithAnInteger* TgtAngCmd;

    G4UIcmdWithAString*        DgrMaterCmd;
    G4UIcmdWithADoubleAndUnit* DgrRadCmd;
    G4UIcmdWithADoubleAndUnit* DgrThickCmd;
    G4UIcmdWithADoubleAndUnit* DgrPosCmd;

    G4UIcmdWithoutParameter*   UpdateCmd;

};

#endif
