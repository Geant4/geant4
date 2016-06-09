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
// $Id: RemSimShieldingDecorator.hh,v 1.6 2005/09/08 06:56:18 guatelli Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
#ifndef RemSimShieldingDecorator_h
#define RemSimShieldingDecorator_h 1

#include "RemSimDecorator.hh"
#include "RemSimVGeometryComponent.hh"
#include "globals.hh"
class G4VPhysicalVolume;
class G4Box;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class RemSimVGeometryComponent;
class RemSimDecorator;
class G4VPhysicalVolume;
class G4VisAttributes;

class RemSimShieldingDecorator: public RemSimDecorator
{
public:
  RemSimShieldingDecorator(RemSimVGeometryComponent*);
  ~RemSimShieldingDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void PrintDetectorParameters();

private:
  void ConstructShielding(G4VPhysicalVolume*);

  G4double shieldingX;
  G4double shieldingY;
  G4double shieldingZ;
  G4double translation;
  RemSimMaterial* pMaterial; 
  G4String shieldingMaterial;
  G4VisAttributes* shieldingVisAtt;
  G4Box* shielding;
  G4LogicalVolume* shieldingLog;
  G4VPhysicalVolume* shieldingPhys;
};
#endif
