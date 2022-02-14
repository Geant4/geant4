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
// GRDetectorConstructionMessenger.hh
//   Headter file of a messenger class that handles geometry configuration.
//
// History
//   September 8th, 2020 : first implementation
//
// ********************************************************************

#ifndef GRDetectorConstructionMessenger_H
#define GRDetectorConstructionMessenger_H 1

#include "G4UImessenger.hh"
#include "globals.hh"

class GRDetectorConstruction;
class G4UIdirectory;
class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

class GRDetectorConstructionMessenger: public G4UImessenger
{
  public:
    GRDetectorConstructionMessenger(GRDetectorConstruction*);
    virtual ~GRDetectorConstructionMessenger();
    virtual void SetNewValue(G4UIcommand*,G4String);
    virtual G4String GetCurrentValue(G4UIcommand*);

  private:
    GRDetectorConstruction* pDC;

    G4UIdirectory* 		geomDir;
    G4UIcmdWithAString* 	selectCmd;
    G4UIcmdWithAnInteger*	listSolidCmd;
    G4UIcmdWithAnInteger*	listLogVolCmd;
    G4UIcmdWithAnInteger*	listPhysVolCmd;
    G4UIcmdWithAnInteger*	listRegionCmd;
    G4UIcommand*                createRegionCmd;
    ////G4UIcommand*                checkOverlapCmd;

    G4UIdirectory*              materialDir;
    G4UIcmdWithAString*         listMatCmd;
    G4UIcmdWithoutParameter*    dumpMatCmd;
    G4UIcmdWithAString*         createMatCmd;
    G4UIcmdWithAString*         getMatCmd;
    G4UIcommand*                setMatCmd;

};

#endif

