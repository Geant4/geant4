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
/// \file optical/LXe/include/LXeDetectorMessenger.hh
/// \brief Definition of the LXeDetectorMessenger class
//
//
#ifndef LXeDetectorMessenger_h
#define LXeDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class LXeDetectorConstruction;

class G4UIcmdWithABool;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAnInteger;
class G4UIcmdWith3VectorAndUnit;
class G4UIcommand;
class G4UIdirectory;

class LXeDetectorMessenger : public G4UImessenger
{
 public:
  LXeDetectorMessenger(LXeDetectorConstruction*);
  ~LXeDetectorMessenger() override;

  void SetNewValue(G4UIcommand*, G4String) override;

 private:
  LXeDetectorConstruction* fLXeDetector = nullptr;
  G4UIdirectory* fDetectorDir = nullptr;
  G4UIdirectory* fVolumesDir = nullptr;
  G4UIcmdWith3VectorAndUnit* fDimensionsCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fHousingThicknessCmd = nullptr;
  G4UIcmdWithADoubleAndUnit* fPmtRadiusCmd = nullptr;
  G4UIcmdWithAnInteger* fNxCmd = nullptr;
  G4UIcmdWithAnInteger* fNyCmd = nullptr;
  G4UIcmdWithAnInteger* fNzCmd = nullptr;
  G4UIcmdWithABool* fSphereCmd = nullptr;
  G4UIcmdWithADouble* fReflectivityCmd = nullptr;
  G4UIcmdWithABool* fWlsCmd = nullptr;
  G4UIcmdWithABool* fLxeCmd = nullptr;
  G4UIcmdWithAnInteger* fNFibersCmd = nullptr;
  G4UIcommand* fDefaultsCmd = nullptr;
  G4UIcmdWithADouble* fMainScintYield = nullptr;
  G4UIcmdWithADouble* fWLSScintYield = nullptr;
  G4UIcmdWithAnInteger* fSaveThresholdCmd = nullptr;
};

#endif
