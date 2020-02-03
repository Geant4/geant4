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
/// \file field/field02/include/F02FieldMessenger.hh
/// \brief Definition of the F02FieldMessenger class
//
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef F02FieldMessenger_h
#define F02FieldMessenger_h 1

#include "G4UImessenger.hh"

class F02ElectricFieldSetup;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithoutParameter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class F02FieldMessenger: public G4UImessenger
{
  public:
    F02FieldMessenger(F02ElectricFieldSetup* );
    virtual ~F02FieldMessenger();

    virtual void SetNewValue(G4UIcommand*, G4String);
 
  private:

    F02ElectricFieldSetup*     fElFieldSetup;

    G4UIdirectory*             fFieldDir;
    G4UIcmdWithAnInteger*      fStepperCmd;
    G4UIcmdWithADoubleAndUnit* fElFieldZCmd;
    G4UIcmdWith3VectorAndUnit* fElFieldCmd;
    G4UIcmdWithADoubleAndUnit* fMinStepCmd;
    G4UIcmdWithoutParameter*   fUpdateCmd;
};

#endif
