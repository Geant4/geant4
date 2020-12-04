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
// GRGeomBiasMessenger.hh
//   Header file of a messenger class that handles the UI commands for
//   geometry imprtance biasing.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRGeomBiasMessenger_h
#define GRGeomBiasMessenger_h 1

class G4ScoringManager;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;

class GRDetectorConstruction;
class GRPhysicsList;

#include "G4UImessenger.hh"

class GRGeomBiasMessenger: public G4UImessenger
{

  public:
    GRGeomBiasMessenger(GRDetectorConstruction*,GRPhysicsList*,G4int verboseLvl = 0);
    ~GRGeomBiasMessenger();
    void SetNewValue(G4UIcommand * command,G4String newValues);
    G4String GetCurrentValue(G4UIcommand * command);

  private:
    GRDetectorConstruction* detector;
    GRPhysicsList* physics;
    G4int verboseLevel;

    G4UIdirectory* biasDir;
    G4UIcommand* geoBiasCmd;
    G4UIcmdWith3VectorAndUnit* geoBiasLocCmd;
    G4UIcmdWithADoubleAndUnit* geoBiasInRadCmd;
    G4UIcmdWith3VectorAndUnit* geoBiasLocTgtCmd;
    G4UIcmdWithAnInteger* geoBiasFucCmd;
    G4UIcmdWithADouble* geoBiasProbCmd;
};

//

#endif
