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
// $Id: ConfigData.hh,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Creation - central configuration data.
//
#ifndef CONFIGDATA_HH
#define CONFIGDATA_HH

#include "globals.hh"
#include <fstream> 

class G4Material;

namespace ConfigData {

  namespace {

    // Geometry configuration with defaults
    G4double targetRadius = 3.63*cm;
    G4double targetDistance = 0.0;
    G4Material* targetMaterial = 0;

    G4double chamberWindowRadius = 2.0*cm;
    G4double chamberWindowDistance = 0.0;
    G4Material* chamberWindowMaterial = 0;

    G4double airGap1Radius = 2.0*cm;
    G4double airGap1Distance = 0.0;
    G4Material* airGap1Material = 0;

    G4double monitorRadius = 2.0*cm;
    G4double monitorDistance = 0.0;
    G4Material* monitorMaterial = 0;

    G4double airGap2Radius = 2.0*cm;
    G4double airGap2Distance = 0.0;
    G4Material* airGap2Material = 0;

    G4double beamWindowRadius = 2.0*cm;
    G4double beamWindowDistance = 0.0;
    G4Material* beamWindowMaterial = 0;

    G4double beamPipeRadius = 5.0*cm;
    G4double beamPipeDistance = 0.0;
    G4Material* beamPipeMaterial = 0;

    G4double beamPipeEndDistance = 0.0;

    // Scorer configuration data with defaults
    G4double scorerMinRadius = 1.0*m;      // Inner scoring sphere radius
    G4double scorerMaxRadius = 1.001*m;    // Outer scoring sphere radius

    G4double scorerMinTheta = 0.0*deg;     // Minimum theta
    G4double scorerMaxTheta = 180.0*deg;   // Maximum theta
    G4double scorerDeltaTheta = 0.5*deg;   // Delta theta

    G4double scorerMinEnergy = 0.0*MeV;    // Minimum energy
    G4double scorerMaxEnergy = 30.0*MeV;   // Maximum energy
    G4double scorerDeltaEnergy = 0.03*MeV; // Delta energy

    // Primary particle configuration
    G4double primaryEnergy = 0.0;

    // MultiFunctionalDetector name
    G4String detectorName = "MultiFunctionalDetector";

    // Output results and configuration dump
    G4String outputDirectory;
    G4String outputId;
    std::ofstream* outFile;

    // Verbose
    G4bool verbose = false;
  }

  void SetTargetRadius(G4double value);
  G4double GetTargetRadius();

  void SetTargetDistance(G4double value);
  G4double GetTargetDistance();

  void SetTargetMaterial(G4Material* value);
  G4Material* GetTargetMaterial();

  void SetChamberWindowRadius(G4double value);
  G4double GetChamberWindowRadius();

  void SetChamberWindowDistance(G4double value);
  G4double GetChamberWindowDistance();

  void SetChamberWindowMaterial(G4Material* value);
  G4Material* GetChamberWindowMaterial();

  void SetAirGap1Radius(G4double value);
  G4double GetAirGap1Radius();

  void SetAirGap1Distance(G4double value);
  G4double GetAirGap1Distance();

  void SetAirGap1Material(G4Material* value);
  G4Material* GetAirGap1Material();

  void SetMonitorRadius(G4double value);
  G4double GetMonitorRadius();

  void SetMonitorDistance(G4double value);
  G4double GetMonitorDistance();

  void SetMonitorMaterial(G4Material* value);
  G4Material* GetMonitorMaterial();

  void SetAirGap2Radius(G4double value);
  G4double GetAirGap2Radius();

  void SetAirGap2Distance(G4double value);
  G4double GetAirGap2Distance();

  void SetAirGap2Material(G4Material* value);
  G4Material* GetAirGap2Material();

  void SetBeamWindowRadius(G4double value);
  G4double GetBeamWindowRadius();

  void SetBeamWindowDistance(G4double value);
  G4double GetBeamWindowDistance();

  void SetBeamWindowMaterial(G4Material* value);
  G4Material* GetBeamWindowMaterial();

  void SetBeamPipeRadius(G4double value);
  G4double GetBeamPipeRadius();

  void SetBeamPipeDistance(G4double value);
  G4double GetBeamPipeDistance();

  void SetBeamPipeMaterial(G4Material* value);
  G4Material* GetBeamPipeMaterial();

  void SetBeamPipeEndDistance(G4double value);
  G4double GetBeamPipeEndDistance();

  ////
  void SetScorerMinRadius(G4double value);
  G4double GetScorerMinRadius();

  void SetScorerMaxRadius(G4double value);
  G4double GetScorerMaxRadius();

  void SetScorerMinTheta(G4double value);
  G4double GetScorerMinTheta();

  void SetScorerMaxTheta(G4double value);
  G4double GetScorerMaxTheta();

  void SetScorerDeltaTheta(G4double value);
  G4double GetScorerDeltaTheta();

  void SetScorerMinEnergy(G4double value);
  G4double GetScorerMinEnergy();

  void SetScorerMaxEnergy(G4double value);
  G4double GetScorerMaxEnergy();

  void SetScorerDeltaEnergy(G4double value);
  G4double GetScorerDeltaEnergy();

  void SetPrimaryEnergy(G4double value);
  G4double GetPrimaryEnergy();

  G4String GetDetectorName(); 

  void SetOutputDirectory(const G4String&);
  G4String GetOutputDirectory();

  void SetOutputId(const G4String&);
  G4String GetOutputId();

  std::ofstream& GetConfigFile();

  void SetVerbose(G4bool);
  G4bool GetVerbose();

  namespace BremSplitting {

    namespace {
      G4bool active = false;
      unsigned factor = 0;

      unsigned nSecondaries = 0;
    }

    G4bool GetActivation();
    unsigned GetFactor();

    void SetActivation(G4bool);
    void SetFactor(unsigned);

    void AddSecondaries(unsigned);
    unsigned GetNSecondaries();
  }
}

#endif
