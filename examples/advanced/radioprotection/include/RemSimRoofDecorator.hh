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
//    ************************************
//    *                                  *
//    *      RemSimRoofDecorator.hh      *
//    *                                  *
//    ************************************
//
// Code developed by: S.Guatelli, guatelli@ge.infn.it
//
//
// $Id: RemSimRoofDecorator.hh,v 1.6 2006-06-29 16:23:13 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef RemSimRoofDecorator_h
#define RemSimRoofDecorator_h 1

#include "RemSimDecorator.hh"
#include "RemSimVGeometryComponent.hh"
#include "globals.hh"
class G4VPhysicalVolume;
class G4Trd;
class G4LogicalVolume;
class G4Material;
class RemSimMaterial;
class RemSimVGeometryComponent;
class RemSimDecorator;
class G4VPhysicalVolume;
class G4VisAttributes;

class RemSimRoofDecorator: public RemSimDecorator
{
public:
  RemSimRoofDecorator(RemSimVGeometryComponent*);
  ~RemSimRoofDecorator();
  void ConstructComponent(G4VPhysicalVolume*);
  void DestroyComponent(); 
  void ChangeThickness(G4double);
  void PrintDetectorParameters();

private:

  void ConstructRoof(G4VPhysicalVolume*);

  RemSimMaterial* pMaterial; 
  G4String shieldingMaterial;
  G4VisAttributes* roofVisAtt;
  G4Trd* roof;
  G4LogicalVolume* roofLog;
  G4VPhysicalVolume* roofPhys;
};
#endif
