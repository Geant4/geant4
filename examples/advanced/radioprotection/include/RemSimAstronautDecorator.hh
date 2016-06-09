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
//    **********************************
//    *                                *
//    *    RemSimAstronautDecorator.hh *
//    *                                *
//    **********************************
//
// $Id: RemSimAstronautDecorator.hh,v 1.7 2005/05/27 14:21:42 guatelli Exp $
//
// Author:Susanna Guatelli, guatelli@ge.infn.it 

#ifndef RemSimAstronautDecorator_h
#define RemSimAstronautDecorator_h 1

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
class RemSimSensitiveDetector;

class RemSimAstronautDecorator: public RemSimDecorator
{
public:
  RemSimAstronautDecorator(RemSimVGeometryComponent*, G4bool);
  ~RemSimAstronautDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void PrintDetectorParameters();

private:
  void ConstructAstronaut(G4VPhysicalVolume*);  
  G4Box* phantom;
  G4LogicalVolume* phantomLog;
  G4VPhysicalVolume* phantomPhys;
  G4bool flag;
  RemSimSensitiveDetector* sensitiveDetector;
};
#endif
