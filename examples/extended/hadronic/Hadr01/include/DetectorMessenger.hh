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
// $Id: DetectorMessenger.hh,v 1.1 2006-06-02 19:00:00 vnivanch Exp $
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

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithABool;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3Vector;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorMessenger: public G4UImessenger
{
public:
  DetectorMessenger(DetectorConstruction* );
  virtual ~DetectorMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  DetectorConstruction* Detector;

  G4UIdirectory*             testemDir;
  G4UIcmdWithAString*        MaterCmd;
  G4UIcmdWithAString*        Mat1Cmd;
  G4UIcmdWithAString*        Mat2Cmd;
  G4UIcmdWithADoubleAndUnit* l1Cmd;
  G4UIcmdWithADoubleAndUnit* l2Cmd;
  G4UIcmdWithADoubleAndUnit* l3Cmd;
  G4UIcmdWithADoubleAndUnit* enCmd;
  G4UIcmdWithADoubleAndUnit* eeCmd;
  G4UIcmdWithADoubleAndUnit* egCmd;
  G4UIcmdWithAnInteger*      NEbinsCmd;
  G4UIcmdWithAnInteger*      NXYbinsCmd;
  G4UIcmdWithAnInteger*      NbinCmd;
  G4UIcmdWithAnInteger*      NumOfAbsCmd;
  G4UIcmdWithAnInteger*      verbCmd;
  G4UIcmdWithoutParameter*   UpdateCmd;
  G4UIcmdWithABool*          ntupCmd;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

