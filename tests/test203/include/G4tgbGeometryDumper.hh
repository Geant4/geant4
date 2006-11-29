#ifndef G4tgbGeometryDumper_HH
#define G4tgbGeometryDumper_HH

#include "globals.hh"
#include <fstream>
#include <map>
#include <vector>
class G4Material;
class G4Element;
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
#include "G4RotationMatrix.hh"

class G4tgbGeometryDumper{

private:
  G4tgbGeometryDumper();

public:
  static G4tgbGeometryDumper* GetInstance();
  ~G4tgbGeometryDumper();

  void DumpGeometry(const G4String& fname );
  G4VPhysicalVolume* GetTopPhysVol();
  void DumpPhysVol( G4VPhysicalVolume* pv );
  void DumpLogVol( G4LogicalVolume* lv );
  void DumpMaterial( G4Material* mat );
  void DumpElement( G4Element* ele);
  void DumpSolid( G4VSolid* solid );
  void DumpBooleanVolume( const G4String& solidType, G4VSolid* so );
  void DumpSolidParams(G4VSolid * so);
  void DumpPolySections(G4int zPlanes, G4double* z, G4double *rmin, G4double *rmax);
  G4String DumpRotationMatrix( G4RotationMatrix* rotm );

private:
  std::vector<G4VPhysicalVolume*> GetPVChildren( G4LogicalVolume* lv );
  G4String GetTGSolidType( G4String& solidtype );
  double MatDeterminant(G4RotationMatrix * ro) ;
  G4double approxTo0( double val );
  G4String itoa(int current);
  G4String AddQuotes( const G4String& str );

  G4bool CheckIfElementExists( const G4String& name, G4Element* );
  G4bool CheckIfMaterialExists( const G4String& name, G4Material* );
  G4bool CheckIfLogVolExists( const G4String& name, G4LogicalVolume* pt );
  G4bool CheckIfSolidExists( const G4String& name, G4VSolid* );
  G4bool CheckIfPhysVolExists( const G4String& name, G4VPhysicalVolume* );
  G4String LookForExistingRotation( const G4RotationMatrix* rotm );
  G4String SupressRefl( G4String name );
  G4String SubstituteRefl( G4String name );

private:
  static G4tgbGeometryDumper* theInstance;

  std::ofstream* theFile;

  std::map<G4String,G4Element*> theElements;
  std::map<G4String,G4Material*> theMaterials;
  std::map<G4String,G4VSolid*> theSolids;
  std::map<G4String,G4LogicalVolume*> theLogVols;
  std::map<G4String,G4VPhysicalVolume*> thePhysVols;
  std::map<G4String,G4RotationMatrix*> theRotMats;

  int theRotationNumber;
};

#endif
