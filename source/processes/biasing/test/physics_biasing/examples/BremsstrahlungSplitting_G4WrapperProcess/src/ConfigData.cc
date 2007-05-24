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
// $Id: ConfigData.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - central configuration data.
//

#include "ConfigData.hh"

namespace ConfigData {

  void SetTargetRadius(G4double value) {targetRadius = value;}
  G4double GetTargetRadius() {return targetRadius;}

  void SetTargetDistance(G4double value) {targetDistance = value;}
  G4double GetTargetDistance() {return targetDistance;}

  void SetTargetMaterial(G4Material* value) {targetMaterial = value;}
  G4Material* GetTargetMaterial() {return targetMaterial;}

  void SetChamberWindowRadius(G4double value) {chamberWindowRadius = value;}
  G4double GetChamberWindowRadius() {return chamberWindowRadius;}

  void SetChamberWindowDistance(G4double value) {chamberWindowDistance = value;}
  G4double GetChamberWindowDistance() {return chamberWindowDistance;}

  void SetChamberWindowMaterial(G4Material* value) {chamberWindowMaterial = value;}
  G4Material* GetChamberWindowMaterial() {return chamberWindowMaterial;}

  void SetAirGap1Radius(G4double value) {airGap1Radius = value;}
  G4double GetAirGap1Radius() {return airGap1Radius;}

  void SetAirGap1Distance(G4double value) {airGap1Distance = value;}
  G4double GetAirGap1Distance() {return airGap1Distance;}

  void SetAirGap1Material(G4Material* value) {airGap1Material = value;}
  G4Material* GetAirGap1Material() {return airGap1Material;}

  void SetMonitorRadius(G4double value) {monitorRadius = value;}
  G4double GetMonitorRadius() {return monitorRadius;}

  void SetMonitorDistance(G4double value) {monitorDistance = value;}
  G4double GetMonitorDistance() {return monitorDistance;}

  void SetMonitorMaterial(G4Material* value) {monitorMaterial = value;}
  G4Material* GetMonitorMaterial() {return monitorMaterial;}

  void SetAirGap2Radius(G4double value) {airGap2Radius = value;}
  G4double GetAirGap2Radius() {return airGap2Radius;}

  void SetAirGap2Distance(G4double value) {airGap2Distance = value;}
  G4double GetAirGap2Distance() {return airGap2Distance;}

  void SetAirGap2Material(G4Material* value) {airGap2Material = value;}
  G4Material* GetAirGap2Material() {return airGap2Material;}

  void SetBeamWindowRadius(G4double value) {beamWindowRadius = value;}
  G4double GetBeamWindowRadius() {return beamWindowRadius;}

  void SetBeamWindowDistance(G4double value) {beamWindowDistance = value;}
  G4double GetBeamWindowDistance() {return beamWindowDistance;}

  void SetBeamWindowMaterial(G4Material* value) {beamWindowMaterial = value;}
  G4Material* GetBeamWindowMaterial() {return beamWindowMaterial;}

  void SetBeamPipeRadius(G4double value) {beamPipeRadius = value;}
  G4double GetBeamPipeRadius() {return beamPipeRadius;}

  void SetBeamPipeDistance(G4double value) {beamPipeDistance = value;}
  G4double GetBeamPipeDistance() {return beamPipeDistance;}

  void SetBeamPipeMaterial(G4Material* value) {beamPipeMaterial = value;}
  G4Material* GetBeamPipeMaterial() {return beamPipeMaterial;}

  void SetBeamPipeEndDistance(G4double value) {beamPipeEndDistance = value;}
  G4double GetBeamPipeEndDistance() {return beamPipeEndDistance;}

  //////
  void SetScorerMinRadius(G4double value) {scorerMinRadius = value;}
  G4double GetScorerMinRadius() {return scorerMinRadius;}

  void SetScorerMaxRadius(G4double value) {scorerMaxRadius = value;}
  G4double GetScorerMaxRadius() {return scorerMaxRadius;}

  void SetScorerMinTheta(G4double value) {scorerMinTheta = value;}
  G4double GetScorerMinTheta() {return scorerMinTheta;}

  void SetScorerMaxTheta(G4double value) {scorerMaxTheta = value;}
  G4double GetScorerMaxTheta() {return scorerMaxTheta;}

  void SetScorerDeltaTheta(G4double value) {scorerDeltaTheta = value;}
  G4double GetScorerDeltaTheta() {return scorerDeltaTheta;}

  void SetScorerMinEnergy(G4double value) {scorerMinEnergy = value;}
  G4double GetScorerMinEnergy() {return scorerMinEnergy;}

  void SetScorerMaxEnergy(G4double value) {scorerMaxEnergy = value;}
  G4double GetScorerMaxEnergy() {return scorerMaxEnergy;}

  void SetScorerDeltaEnergy(G4double value) {scorerDeltaEnergy = value;}
  G4double GetScorerDeltaEnergy() {return scorerDeltaEnergy;}

  // 
  void SetPrimaryEnergy(G4double value) {primaryEnergy = value;}
  G4double GetPrimaryEnergy() {return primaryEnergy;}

  G4String GetDetectorName() {return detectorName;}
  
  void SetOutputDirectory(const G4String& dir) {outputDirectory = dir;}
  G4String GetOutputDirectory() {return outputDirectory;}

  void SetOutputId(const G4String& id) {outputId = id;}
  G4String GetOutputId() {return outputId;}

  void SetVerbose(G4bool value) {verbose = value;}
  G4bool GetVerbose() {return verbose;}

  std::ofstream& GetConfigFile() 
  {
    if (0 != outFile) return *outFile;

    outFile = new std::ofstream(G4String(outputDirectory+"/"+outputId+".Config").data());

    return *outFile;
  }

  namespace BremSplitting {
    G4bool GetActivation() {return active;}
    unsigned GetFactor() {return factor;}

    void SetActivation(G4bool value) {active = value;}
    void SetFactor(unsigned value) {factor = value;}

    void AddSecondaries(unsigned value) {nSecondaries += value;}
    unsigned GetNSecondaries() {return nSecondaries;}
  }
}

