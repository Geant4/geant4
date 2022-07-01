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
// G4ParticleGunMessenger
//
// Class description:
//
// This is a concrete class of G4UImessenger which handles commands for
// G4ParticleGun.

// Author: Makoto Asai, 1997
// --------------------------------------------------------------------
#ifndef G4ParticleGunMessenger_hh
#define G4ParticleGunMessenger_hh 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4ParticleGun;
class G4ParticleTable;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithAnInteger;

class G4ParticleGunMessenger : public G4UImessenger
{
  public:

    explicit G4ParticleGunMessenger(G4ParticleGun* fPtclGun);
    ~G4ParticleGunMessenger() override;
    
    void SetNewValue(G4UIcommand* command, G4String newValues) override;
    G4String GetCurrentValue(G4UIcommand* command) override;

  private:

    void IonCommand(const G4String& newValues);
    void IonLevelCommand(const G4String& newValues);

  private:

    G4ParticleGun* fParticleGun = nullptr; // Not owned, cannot be null post construction
    G4ParticleTable* particleTable = nullptr;

    // Commands
    //
    G4UIdirectory*              gunDirectory;
    G4UIcmdWithoutParameter*    listCmd;
    G4UIcmdWithAString*         particleCmd;
    G4UIcmdWith3Vector*         directionCmd;
    G4UIcmdWithADoubleAndUnit*  energyCmd;
    G4UIcmdWithADoubleAndUnit*  momAmpCmd;
    G4UIcmdWith3VectorAndUnit*  momCmd;
    G4UIcmdWith3VectorAndUnit*  positionCmd;
    G4UIcmdWithADoubleAndUnit*  timeCmd;
    G4UIcmdWith3Vector*         polCmd;
    G4UIcmdWithAnInteger*       numberCmd;
    G4UIcommand*                ionCmd;
    G4UIcommand*                ionLvlCmd;

    // For ion shooting
    //
    G4bool   fShootIon = false; 
    G4int    fAtomicNumber = 0;
    G4int    fAtomicMass = 0;
    G4int    fIonCharge = 0;
    G4double fIonExciteEnergy = 0.0;
    char     fIonFloatingLevelBase = '\0';
    G4int    fIonEnergyLevel = 0;
};

#endif
