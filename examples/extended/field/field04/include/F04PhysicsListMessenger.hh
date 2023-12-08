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
/// \file field/field04/include/F04PhysicsListMessenger.hh
/// \brief Definition of the F04PhysicsListMessenger class
//

#ifndef F04PhysicsListMessenger_h
#define F04PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class F04PhysicsList;

class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

/// Provide control of the physics list and cut parameters

class F04PhysicsListMessenger : public G4UImessenger
{
public:

    F04PhysicsListMessenger(F04PhysicsList* );
    ~F04PhysicsListMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

private:

    F04PhysicsList* fPhysicsList = nullptr;

    G4UIdirectory* fDirectory = nullptr;
    G4UIdirectory* fDecayDirectory = nullptr;

    G4UIcmdWithADoubleAndUnit* fStepMaxCMD = nullptr;
/*
    G4UIcmdWithAString*        fRemovePhysicsCMD = nullptr;
    G4UIcmdWithoutParameter*   fClearPhysicsCMD = nullptr;
*/
    G4UIcmdWithoutParameter* fPienuCMD = nullptr;
    G4UIcmdWithoutParameter* fPimunuCMD = nullptr;
};
#endif
