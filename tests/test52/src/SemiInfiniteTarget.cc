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

#include "SemiInfiniteTarget.hh"
#include "TargetGeometryManager.hh"
#include "TargetLayerSD.hh"
#include "Materials.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"


SemiInfiniteTarget::SemiInfiniteTarget(G4String layerName,
                                       TargetGeometryManager* geomManager,
                                       G4VPhysicalVolume* world,
                                       G4String material,
                                       G4double thickn,
                                       G4double rad,
                                       G4double maxStep) :
    TargetComponent(geomManager, layerName, thickn, material),
    worldVolPhys(world), 
    radius(rad),
    maxStepSize(maxStep),
    semiInfTargetSolid(0),
    semiInfTargetLogic(0),
    semiInfTargetPhys(0) {

  updateManager -> Attach(this);  

  semiInfTargetSolid = new G4Tubs(Name(), 0.0 * mm, radius, Thickness() * 0.5,
				  0.0 * deg, 360.0 * deg); 
      
  semiInfTargetLogic = new G4LogicalVolume(semiInfTargetSolid,
                                           Material(),
                                           Name());

  semiInfTargetPhys = new G4PVPlacement(0, 
                         G4ThreeVector(0.0 * cm, 0.0 * cm, Thickness() * 0.5),
                         Name(), semiInfTargetLogic, worldVolPhys, false, 0);

  TargetLayerSD* semiInfTargetSD = new TargetLayerSD(Name()); 

  semiInfTargetLogic -> SetSensitiveDetector(semiInfTargetSD);

  G4SDManager* SDManager = G4SDManager::GetSDMpointer();
  SDManager -> AddNewDetector(semiInfTargetSD);

  semiInfTargetUserLimits = new G4UserLimits(maxStepSize);
  semiInfTargetLogic -> SetUserLimits(semiInfTargetUserLimits);

  G4VisAttributes* semiInfTargetVisAtt = 
                      new G4VisAttributes(G4Colour::Yellow());
  semiInfTargetVisAtt -> SetVisibility(true);
  semiInfTargetLogic -> SetVisAttributes(semiInfTargetVisAtt);

  GeometryUpdate(updateManager);
}


SemiInfiniteTarget::~SemiInfiniteTarget() {

}


void SemiInfiniteTarget::GeometryUpdate(TargetGeometryManager* geomManager) {

  if(geomManager == updateManager) {

     G4double newRadius = updateManager -> GetRadius();
     if(newRadius > 0.0 * cm) {
	radius = newRadius;
        semiInfTargetSolid -> SetOuterRadius(radius);
     } 

     G4double newMaxStepSize = updateManager -> GetMaxStepSize();
     if(newMaxStepSize > 0.0 * mm) {
        maxStepSize = newMaxStepSize;
        semiInfTargetUserLimits -> SetMaxAllowedStep(maxStepSize);
     }

     Thickness(updateManager -> GetThickness(this));
     semiInfTargetSolid -> SetZHalfLength(Thickness() * 0.5);

     Material(updateManager -> GetMaterial(this));
     semiInfTargetLogic -> SetMaterial(Material());     

     G4double position = updateManager -> GetPosition(this);
     semiInfTargetPhys -> SetTranslation(
         G4ThreeVector(0.0 * cm, 0.0 * cm, position));
  }
}


