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
// $Id: TiaraGeometry.hh,v 1.5 2006/06/29 15:43:46 gunter Exp $
// GEANT4 tag $Name: geant4-09-00 $
//
// ----------------------------------------------------------------------
//
// Class TiaraGeometry
//

#ifndef TiaraGeometry_hh
#define TiaraGeometry_hh TiaraGeometry_hh

#include <set>

#include "TiaraDimensions.hh"
#include "G4VUserDetectorConstruction.hh"
#include "TiaraVComponent.hh"

typedef std::set<TiaraVComponent*> TiaraComponents;

class G4Material;
class G4LogicalVolume;
class TiaraMaterials;

class TiaraGeometry : public G4VUserDetectorConstruction {
public:
  TiaraGeometry(TiaraMaterials &matfac);
  ~TiaraGeometry();

  void BuildGeometry(const TiaraDimensions &td);
  
  virtual G4VPhysicalVolume* Construct();
  
  G4Material *GetWorldMaterial();

  G4VPhysicalVolume *PlaceExpComponent(const G4ThreeVector &pos, 
				       G4LogicalVolume *logVol,
				       const G4String &physName);  


  G4LogicalVolume *BuildCollimator(G4double width,
				   const G4String &outerMatName,
				   const G4String &innerMatName);


  G4LogicalVolume *BuildShield(G4double width,
			       const G4String &matName);

  G4VPhysicalVolume *AddPhysicalDetector(G4double xDist,
					 const G4String &physName);

  G4VPhysicalVolume *AddPhysicalRingDetector(G4double xDist,
					     const G4String &physName);
  G4VPhysicalVolume *AddDetectorSlab( const G4String &physName);
  G4VPhysicalVolume *AddSourceDetector();
  void CreateComponents();

  
private:
  void ConstHall();
  void PlaceComponents();

  G4LogicalVolume *fLogicalWorld;
  G4VPhysicalVolume *fWorldVolume;
  TiaraMaterials &fMaterials;
  TiaraComponents fTiaraComponents;
  G4double fShieldWidth;
  G4double fColimatorWidth;
  G4LogicalVolume *fColliPipeLogVol;
  TiaraDimensions fTiaraDimensions;
};


#endif
