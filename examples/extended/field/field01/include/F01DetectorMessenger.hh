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
/// \file field/field01/include/F01DetectorMessenger.hh
/// \brief Definition of the F01DetectorMessenger class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F01DetectorMessenger_h
#define F01DetectorMessenger_h 1

#include "G4UImessenger.hh"

class F01DetectorConstruction;

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F01DetectorMessenger: public G4UImessenger
{
  public:

    F01DetectorMessenger(F01DetectorConstruction* );
    virtual ~F01DetectorMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);

  private:

    F01DetectorConstruction*   fDetector;

    G4UIdirectory*             fDetDir;

    G4UIcmdWithAString*        fAbsMaterCmd;
    G4UIcmdWithADoubleAndUnit* fAbsThickCmd;
    G4UIcmdWithADoubleAndUnit* fAbsRadCmd;

    G4UIcmdWithADoubleAndUnit* fAbsZposCmd;

    G4UIcmdWithAString*        fWorldMaterCmd;
    G4UIcmdWithADoubleAndUnit* fWorldZCmd;
    G4UIcmdWithADoubleAndUnit* fWorldRCmd;

};

#endif
