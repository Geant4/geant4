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
#ifndef DETECTORCONSTRUCTION_HH
#define DETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectorMessenger;
class TargetGeometryManager;
class AnalysisBuilder;
class G4VPhysicalVolume;


class DetectorConstruction : public G4VUserDetectorConstruction {

 public:
   DetectorConstruction();
   ~DetectorConstruction();

   void CreateFrontLayer(G4String layerName);
   void SetLayerRadius(G4double rad);
   void SetLayerThickness(G4double thickn);
   void SetLayerMaterial(G4String mat);
   void SetLayerMaxStepSize(G4double max);
   void CreateCalorimeter(G4double zPosition);
   void SetCalorimeterThickness(G4double thickn);

 private:
   G4VPhysicalVolume* Construct();    

   DetectorMessenger* messenger;
   TargetGeometryManager* updateManager;
 
   G4VPhysicalVolume* worldVolPhys;
   G4String targetRegionName;
   G4double calThickness;
};

#endif // DETECTORCONSTRUCTION_HH
