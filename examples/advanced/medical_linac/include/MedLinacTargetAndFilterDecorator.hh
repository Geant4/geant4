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
// $Id: MedLinacTargetAndFilterDecorator.hh,v 1.2 2005/07/03 23:27:37 mpiergen Exp $
//
// Code developed by: M. Piergentili
//
//


#ifndef MedLinacTargetAndFilterDecorator_h
#define MedLinacTargetAndFilterDecorator_h 1

#include "MedLinacDecorator.hh"
#include "MedLinacVGeometryComponent.hh"
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class G4VisAttributes;
class MedLinacVGeometryComponent;
class MedLinacDecorator;

class MedLinacTargetAndFilterDecorator: public MedLinacDecorator
{
public:
  MedLinacTargetAndFilterDecorator(MedLinacVGeometryComponent*);
  ~MedLinacTargetAndFilterDecorator();
  void ConstructComponent(G4VPhysicalVolume*,G4VPhysicalVolume*);
  void DestroyComponent(); 

private:
  void ConstructTargetAndFilter(G4VPhysicalVolume*,G4VPhysicalVolume*);
  G4LogicalVolume* targetA_log;
  G4LogicalVolume* targetB_log;
  G4LogicalVolume* layer1_log;
  G4LogicalVolume* layer2_log;
  G4LogicalVolume* layer3_log;
  G4LogicalVolume* layer4_log;
  G4LogicalVolume* layer5_log;
  G4LogicalVolume* layer6_log;
  G4LogicalVolume* layer7_log;
  G4LogicalVolume* layer8_log;
  G4LogicalVolume* layer9_log;
  G4LogicalVolume* layer10_log;
  G4LogicalVolume* layer11_log;
  G4LogicalVolume* layer12_log;
  G4LogicalVolume* layer13_log;
  G4LogicalVolume* layer14_log;
  G4LogicalVolume* layer15_log;
  G4LogicalVolume* layer16_log;
  G4LogicalVolume* layer17_log;
  G4LogicalVolume* layer18_log;
  G4LogicalVolume* layer19_log;
  G4LogicalVolume* layer20A_log;
  G4LogicalVolume* cone20_log;
  G4LogicalVolume* layer20_log;
  G4LogicalVolume* layer21_log;

  G4VPhysicalVolume* targetA_phys;
  G4VPhysicalVolume* targetB_phys;
  G4VPhysicalVolume* layer1_phys;
  G4VPhysicalVolume* layer2_phys;
  G4VPhysicalVolume* layer3_phys;
  G4VPhysicalVolume* layer4_phys;
  G4VPhysicalVolume* layer5_phys;
  G4VPhysicalVolume* layer6_phys;
  G4VPhysicalVolume* layer7_phys;
  G4VPhysicalVolume* layer8_phys;
  G4VPhysicalVolume* layer9_phys;
  G4VPhysicalVolume* layer10_phys;
  G4VPhysicalVolume* layer11_phys;
  G4VPhysicalVolume* layer12_phys;
  G4VPhysicalVolume* layer13_phys;
  G4VPhysicalVolume* layer14_phys;
  G4VPhysicalVolume* layer15_phys;
  G4VPhysicalVolume* layer16_phys;
  G4VPhysicalVolume* layer17_phys;
  G4VPhysicalVolume* layer18_phys;
  G4VPhysicalVolume* layer19_phys;
  G4VPhysicalVolume* layer20_phys;
  G4VPhysicalVolume* layer21_phys;


  G4VisAttributes* simpleTungstenSVisAtt;
  G4VisAttributes* simpleCopperSVisAtt;
  G4VisAttributes* simpleWorldVisAtt;
  //MedLinacVGeometryComponent* component;
};
#endif
