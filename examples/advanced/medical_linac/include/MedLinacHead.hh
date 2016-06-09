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
// $Id: MedLinacHead.hh,v 1.3 2006/06/29 16:03:43 gunter Exp $
//
// Code developed by: M. Piergentili
//
//
#ifndef MedLinacHead_h
#define MedLinacHead_h 1

class MedLinacVGeometryComponent;
class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Region;
class G4Material;
class G4VisAttributes;
class MedLinacHead: public MedLinacVGeometryComponent
{
public:
  MedLinacHead();
  ~MedLinacHead();
  void ConstructComponent(G4VPhysicalVolume*,G4VPhysicalVolume*);
  void DestroyComponent(); 
  
  private:
  G4LogicalVolume* windowUp_log;
  G4LogicalVolume* UpperCollimator_log;
  G4LogicalVolume* collim_log;
  G4LogicalVolume* tracker_log;
  G4LogicalVolume* CylMinusCone_log;
  G4LogicalVolume* windowLow_log;
  G4LogicalVolume* Window_log;
  G4LogicalVolume* SignalPlate_log;
  G4LogicalVolume* Mirror_log;
  G4LogicalVolume* reticle_log;

  G4VPhysicalVolume* windowUp_phys;
  G4VPhysicalVolume* UpperCollimator_phys;
  G4VPhysicalVolume* CylMinusCone_phys;
  G4VPhysicalVolume* windowLow_phys;
  G4VPhysicalVolume* Window1_phys;
  G4VPhysicalVolume* SignalPlate1_phys;
  G4VPhysicalVolume* SignalPlate2_phys;
  G4VPhysicalVolume* Window2_phys;
  G4VPhysicalVolume* SignalPlate3_phys;
  G4VPhysicalVolume* SignalPlate4_phys;
  G4VPhysicalVolume* Window3_phys;
  G4VPhysicalVolume* Mirror_phys;
  G4VPhysicalVolume* reticle_phys;

  G4Region* aLowerCollRegion;
  G4Region* aUpperCollRegion;


  G4VisAttributes* simpleTungstenWVisAtt;
  G4VisAttributes* simpleTungstenSVisAtt;
  G4VisAttributes* simpleCopperSVisAtt;
  G4VisAttributes* simpleMylarVisAtt;
  G4VisAttributes* simpleKaptonVisAtt;
  G4VisAttributes* simpleWorldVisAtt;
};
#endif
