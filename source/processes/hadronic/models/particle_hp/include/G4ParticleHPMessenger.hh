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
#ifndef G4ParticleHPMessenger_h
#define G4ParticleHPMessenger_h

#include "G4UImessenger.hh"
#include "globals.hh"

class G4ParticleHPManager;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;

class G4ParticleHPMessenger : public G4UImessenger
{
  public:
    G4ParticleHPMessenger(G4ParticleHPManager*);
    ~G4ParticleHPMessenger() override;

    void SetNewValue(G4UIcommand*, G4String) override;

  private:
    G4ParticleHPManager* manager;

    G4UIdirectory* ParticleHPDir;
    G4UIcmdWithAString* PhotoEvaCmd;
    G4UIcmdWithAString* SkipMissingCmd;
    G4UIcmdWithAString* NeglectDopplerCmd;
    G4UIcmdWithAString* DoNotAdjustFSCmd;
    G4UIcmdWithAString* ProduceFissionFragementCmd;
    G4UIcmdWithAString* WendtFissionModelCmd;
    G4UIcmdWithAString* NRESP71Cmd;
    G4UIcmdWithAnInteger* VerboseCmd;
    G4UIcmdWithABool* UseDBRCCmd;
    G4UIcmdWithADoubleAndUnit* MaxEnergySVTCmd;
    G4UIcmdWithADouble* MinADBRCCmd;
    G4UIcmdWithADoubleAndUnit* MinEnergyDBRCCmd;
    G4UIcmdWithADoubleAndUnit* MaxEnergyDBRCCmd;
};

#endif
