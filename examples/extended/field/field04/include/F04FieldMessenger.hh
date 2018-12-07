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
/// \file field/field04/include/F04FieldMessenger.hh
/// \brief Definition of the F04FieldMessenger class
//

#ifndef F04FieldMessenger_h
#define F04FieldMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class F04GlobalField;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class F04DetectorConstruction;

class F04FieldMessenger: public G4UImessenger
{
  public:
    F04FieldMessenger(F04GlobalField*, F04DetectorConstruction* );
    virtual ~F04FieldMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:

    F04GlobalField*            fGlobalField;
 
    G4UIdirectory*             fDetDir;

    G4UIcmdWithADoubleAndUnit* fCaptureB1Cmd;
    G4UIcmdWithADoubleAndUnit* fCaptureB2Cmd;
    G4UIcmdWithADoubleAndUnit* fTransferBCmd;

    G4UIcmdWithAnInteger*      fStepperCMD;
    G4UIcmdWithADoubleAndUnit* fMinStepCMD;
    G4UIcmdWithADoubleAndUnit* fDeltaChordCMD;
    G4UIcmdWithADoubleAndUnit* fDeltaOneStepCMD;
    G4UIcmdWithADoubleAndUnit* fDeltaIntersectionCMD;
    G4UIcmdWithADoubleAndUnit* fEpsMinCMD;
    G4UIcmdWithADoubleAndUnit* fEpsMaxCMD;

    F04DetectorConstruction*   fDetector;
};

#endif
