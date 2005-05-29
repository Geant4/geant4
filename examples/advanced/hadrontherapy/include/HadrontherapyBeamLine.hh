//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: HadrontherapyDetectorConstruction.hh,v 3.0, September 2004
// --------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// --------------------------------------------------------------
// Code developed by:
//
// G. Candiano, G.A.P. Cirrone, F. Di Rosa, G. Russo
// Laboratori Nazionali del Sud - INFN, Catania, Italy
//
// --------------------------------------------------------------
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
  void HadrontherapyBeamScatteringFoils();
  void HadrontherapyBeamCollimators();
  void HadrontherapyBeamMonitoring();
  void HadrontherapyBeamNozzle();
  void HadrontherapyBeamFinalCollimator();
  void setRangeShifterXPos(G4double);
  void setRangeShifterX(G4double);
  void SetFirstScatteringFoil(G4double);
  void SetSecondScatteringFoil(G4double);
  void SetOuterRadiusStopper(G4double);
  void SetInnerRadiusFinalCollimator(G4double);
  void SetRSMaterial(G4String);

private:
  G4VPhysicalVolume* physiBeamLineSupport; 
  G4VPhysicalVolume* physiBeamLineCover; 
  G4VPhysicalVolume* physiBeamLineCover2;

  G4Box* FirstScatteringFoil;
  G4VPhysicalVolume* physiFirstScatteringFoil;
  G4VPhysicalVolume* physiVacuumZone;
  G4VPhysicalVolume* physiKaptonWindow;

  G4Tubs* solidStopper;

  G4VPhysicalVolume* physiStopper; 
  G4Box* SecondScatteringFoil;  
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
   
  G4double FirstScatteringFoil_x;
  G4double outerRadiusStopper;
  G4double SecondScatteringFoil_x;
  G4double RangeShifterBox_x;
  G4double RangeShifterBoxPosition_x; 
  G4double innerRadiusFinalCollimator;

  G4double RangeShifterBox_y; 
  G4double RangeShifterBox_z;
  G4double RangeShifterBoxPosition_y;
  G4double RangeShifterBoxPosition_z;

  G4VPhysicalVolume* mother;
  HadrontherapyMaterial* material;
  G4Material* RSMat;
};
#endif
