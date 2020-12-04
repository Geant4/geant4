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
//  Gorad (Geant4 Open-source Radiation Analysis and Design)
//
//  Author : Makoto Asai (SLAC National Accelerator Laboratory)
//
//  Development of Gorad is funded by NASA Johnson Space Center (JSC)
//  under the contract NNJ15HK11B.
//
// ********************************************************************
//
// GRPhysicsListMessenger.hh
//   Header file of a messenger class that handles Gorad physics list
//   options.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRPhysicsListMessenger_H
#define GRPhysicsListMessenger_H 1

#include "G4UImessenger.hh"
#include "globals.hh"

class GRPhysicsList;
class G4UIcommand;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;

class GRPhysicsListMessenger: public G4UImessenger
{
  public:
    GRPhysicsListMessenger(GRPhysicsList*);
    virtual ~GRPhysicsListMessenger();
    virtual void SetNewValue(G4UIcommand*,G4String);
    virtual G4String GetCurrentValue(G4UIcommand*);

  private:
    GRPhysicsList* pPL;
    G4UIdirectory* physDir;
    G4UIcmdWithAString* selectEMCmd;
    G4UIcmdWithAString* selectHadCmd;
    G4UIcmdWithoutParameter*   addHPCmd;
    G4UIcmdWithoutParameter*   addRDMCmd;
    G4UIcmdWithoutParameter*   addRMCCmd;
    G4UIcmdWithoutParameter*   addOpticalCmd;
    G4UIcmdWithAString* addStepLimitCmd;

    G4UIdirectory* physLimitDir;
    G4UIcmdWithADoubleAndUnit* setStepLimitCmd;
    G4UIcommand*        setRegionStepLimitCmd;

    G4UIdirectory* physCutDir;
    G4UIcmdWithADoubleAndUnit* setCutCmd;
    G4UIcommand*        setCutParticleCmd;
    G4UIcommand*        setCutRegionCmd;
    G4UIcommand*        setCutRegionParticleCmd;

};

#endif

