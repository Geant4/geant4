// $Id: TiaraGeometry.hh,v 1.3 2003-06-18 16:40:23 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
