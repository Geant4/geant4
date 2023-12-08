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
/// \file field/field04/include/F04DetectorMessenger.hh
/// \brief Definition of the F04DetectorMessenger class
//

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
    ~F04DetectorMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:

    F04DetectorConstruction*   fDetector = nullptr;

    G4UIdirectory*          fDetDir = nullptr;

    G4UIcmdWithAString*        fWorldMaterCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fWorldRCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fWorldZCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* fCaptureRCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fCaptureZCmd = nullptr;

    G4UIcmdWithADoubleAndUnit* fTransferRCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fTransferZCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fTransferPCmd = nullptr;

    G4UIcmdWithAString*        fTgtMaterCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fTgtRadCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fTgtThickCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fTgtPosCmd = nullptr;
    G4UIcmdWithAnInteger*      fTgtAngCmd = nullptr;

    G4UIcmdWithAString*        fDgrMaterCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fDgrRadCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fDgrThickCmd = nullptr;
    G4UIcmdWithADoubleAndUnit* fDgrPosCmd = nullptr;
};

#endif
