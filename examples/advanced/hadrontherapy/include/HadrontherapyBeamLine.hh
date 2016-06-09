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
// ----------------------------------------------------------------------------
// $Id: HadrontherapyBeamLine.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#ifndef HadrontherapyBeamLine_H
#define HadrontherapyBeamLine_H 1

class G4VPhysicalVolume;
class HadrontherapyMaterial;

class HadrontherapyBeamLine 
{
public:
  HadrontherapyBeamLine(G4VPhysicalVolume*);
  ~HadrontherapyBeamLine();

  void HadrontherapyBeamLineSupport();
  // Definition of the beam line support

  void HadrontherapyBeamScatteringFoils();
  // Definition of the first scattering foil, 
  // of the Kapton window, of the stopper 

  void HadrontherapyBeamCollimators();
  // Definition of the first collimator, of the range shifter, 
  // of the second collimator, of the first and second 
  // collimator modulators
 
  void HadrontherapyBeamMonitoring();
  // Definition of three monitor chambers

  void HadrontherapyBeamNozzle();
  // Definition of the beam noozle

  void HadrontherapyBeamFinalCollimator();
  // Definition of the final collimator

  // The following methods allow to change parameters
  // of some beam line components

  void SetRangeShifterXPosition(G4double value);
  // This method allows to move the Range Shifter along
  // the X axis

  void SetRangeShifterXSize(G4double halfSize);
  // This method allows to change the size of the range shifter along
  // the X axis
 
  void SetFirstScatteringFoilXSize(G4double);
  // This method allows to change the size of the first scattering foil
  // along the X axis

  void SetSecondScatteringFoilXSize(G4double);
  // This method allows to change the size of the second scattering foil
  // along the X axis 
  
  void SetOuterRadiusStopper(G4double);
  // This method allows to change the size of the outer radius of the stopper
 
  void SetInnerRadiusFinalCollimator(G4double);
  // This method allows to change the size of the inner radius of the 
  // final collimator
  
  void SetRSMaterial(G4String);
  // This method allows to change the material 
  // of the range shifter

private:
  G4VPhysicalVolume* physiBeamLineSupport; 
  G4VPhysicalVolume* physiBeamLineCover; 
  G4VPhysicalVolume* physiBeamLineCover2;
  G4Box* firstScatteringFoil;
  G4VPhysicalVolume* physiFirstScatteringFoil;
  G4VPhysicalVolume* physiKaptonWindow;
  G4Tubs* solidStopper;
  G4VPhysicalVolume* physiStopper; 
  G4Box* secondScatteringFoil;  
  G4VPhysicalVolume* physiSecondScatteringFoil;  
  G4VPhysicalVolume* physiFirstCollimator;  
  G4VPhysicalVolume* physiHoleFirstCollimator;
  G4Box* solidRangeShifterBox;
  G4LogicalVolume* logicRangeShifterBox;
  G4VPhysicalVolume* physiRangeShifterBox;
  G4VPhysicalVolume* physiSecondCollimator;
  G4VPhysicalVolume* physiHoleSecondCollimator; 
  G4VPhysicalVolume* physiFirstCollimatorModulatorBox;
  G4VPhysicalVolume* physiHoleFirstCollimatorModulatorBox; 
  G4VPhysicalVolume* physiSecondCollimatorModulatorBox;
  G4VPhysicalVolume* physiHoleSecondCollimatorModulatorBox;
  G4VPhysicalVolume* physiFirstMonitorLayer1;
  G4VPhysicalVolume* physiFirstMonitorLayer2;
  G4VPhysicalVolume* physiFirstMonitorLayer3;
  G4VPhysicalVolume* physiFirstMonitorLayer4;
  G4VPhysicalVolume* physiSecondMonitorLayer1;
  G4VPhysicalVolume* physiSecondMonitorLayer2;
  G4VPhysicalVolume* physiSecondMonitorLayer3;
  G4VPhysicalVolume* physiSecondMonitorLayer4;
  G4VPhysicalVolume* physiThirdMonitorLayer1;
  G4VPhysicalVolume* physiThirdMonitorLayer2;
  G4VPhysicalVolume* physiThirdMonitorLayer3;
  G4VPhysicalVolume* physiThirdMonitorLayer4;
  G4VPhysicalVolume* physiNozzleSupport;
  G4VPhysicalVolume* physiHoleNozzleSupport; 
  G4VPhysicalVolume* physiSecondHoleNozzleSupport;
  G4Tubs* solidFinalCollimator; 
  G4VPhysicalVolume* physiFinalCollimator; 
  G4double firstScatteringFoilXSize;
  G4double outerRadiusStopper;
  G4double secondScatteringFoilXSize;
  G4double rangeShifterXSize;
  G4double rangeShifterXPosition;
  G4double innerRadiusFinalCollimator;
  G4VPhysicalVolume* mother;
  HadrontherapyMaterial* material;
  G4Material* RSMat;
};
#endif
