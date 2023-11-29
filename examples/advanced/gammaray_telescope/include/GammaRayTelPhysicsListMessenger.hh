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
// Authors: Susanna Guatelli, susanna@uow.edu.au,
// Authors: Jeremy Davis, jad028@uowmail.edu.au
//

// Code based on the hadrontherapy && radioprotection advanced example

#ifndef GammaRayTelPhysicsListMessenger_h
#define GammaRayTelPhysicsListMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

#include <memory>

class GammaRayTelPhysicsList;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;

class GammaRayTelPhysicsListMessenger: public G4UImessenger {
public:
    explicit GammaRayTelPhysicsListMessenger(GammaRayTelPhysicsList*);

    ~GammaRayTelPhysicsListMessenger() override;

    void SetNewValue(G4UIcommand *command, G4String newValue) override;

private:
    GammaRayTelPhysicsList *pPhysicsList;

    G4UIdirectory *physDir;

    G4UIcmdWithADoubleAndUnit *gammaCutCmd;

    G4UIcmdWithADoubleAndUnit *electronCutCmd;

    G4UIcmdWithADoubleAndUnit *protonCutCmd;

    G4UIcmdWithADoubleAndUnit *allCutCmd;

    G4UIcmdWithAString *pListCmd;

    G4UIcmdWithAString *packageListCmd;
};
#endif
