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
// $Id: TiaraGeometry.hh,v 1.4 2003/06/25 09:12:42 gunter Exp $
// GEANT4 tag $Name: geant4-08-00 $
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
