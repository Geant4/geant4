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
//
// $Id: MedLinacHead.hh,v 1.1 2004/05/14 18:25:39 mpiergen Exp $
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


  G4VisAttributes* simpleLeadWVisAtt;
  G4VisAttributes* simpleLeadSVisAtt;
  G4VisAttributes* simpleAlVisAtt;
  G4VisAttributes* simpleCopperSVisAtt;
  G4VisAttributes* simpleMylarVisAtt;
  G4VisAttributes* simpleKaptonVisAtt;
  G4VisAttributes* simpleWorldVisAtt;
};
#endif
