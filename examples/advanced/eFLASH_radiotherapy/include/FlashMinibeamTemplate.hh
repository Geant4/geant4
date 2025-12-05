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
//

#ifndef FlashMinibeamTemplate_H
#define FlashMinibeamTemplate_H 1

#include "FlashDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4Region;
class G4UserLimits;
class FlashMinibeamTemplate {
public:
  FlashMinibeamTemplate(G4VPhysicalVolume *, G4double, G4double);
  ~FlashMinibeamTemplate();
  void ConstructCollimator();
  G4double hight;
private:
G4double maxStep;
G4Region *TemplateRegion;
G4UserLimits* fStepLimit;
  G4VPhysicalVolume *motherPhys;
  G4double X_center;
   G4double y_center;
    G4double z_center;
     G4double d_between_holes; 
  G4double hole_side ;
  G4double radius;
  G4double field_side;
  G4Material * TEFLON;
    G4Material * TUNGSTEN;
  void SetMaterial();
  void ConstructColl_full();
  void ConstructColl_template();
  void Construct_hole();
void ConstructColl_template_planar();
  G4VisAttributes *blue;
    G4VisAttributes *red;

  // G4LogicalVolume * LeafLV;
};
#endif
