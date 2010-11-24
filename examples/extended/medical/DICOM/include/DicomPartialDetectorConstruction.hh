#ifndef DicomPartialDetectorConstruction_h
#define DicomPartialDetectorConstruction_h 1

#include "globals.hh"
#include <map>
#include "DicomDetectorConstruction.hh"

class G4PartialPhantomParameterisation;
class G4LogicalVolume;
class G4Material;


struct matInfo 
{
  G4double sumdens;
  G4int nvoxels;
  G4int id;
};

class DicomPartialDetectorConstruction : public DicomDetectorConstruction
{
public:

  DicomPartialDetectorConstruction();
  ~DicomPartialDetectorConstruction();

  G4VPhysicalVolume* Construct();

private:
  virtual void ReadPhantomData();
  void ReadPhantomDataFile(const G4String& fname);
  void ConstructPhantomContainer();
  void ConstructPhantom();

  void ReadVoxelDensitiesPartial( std::ifstream& fin, std::map< G4int, std::map< G4int, G4int > > ifxmin, std::map< G4int, std::map< G4int, G4int > > ifxmax );

  std::pair<G4double,G4double> ReadVoxelDim( G4int nVoxel, std::ifstream& fin ); 
  void ReadVoxelDensitiesPartial( std::ifstream& fin );

  void SetScorer(G4LogicalVolume* voxel_logic);

  G4Material* BuildMaterialChangingDensity( const G4Material* origMate, float density, G4String mateName );


private:
  G4PartialPhantomParameterisation* thePartialPhantomParam;

  std::multimap<G4int,G4int> theFilledIDs;
  std::map< G4int, std::map< G4int, G4int > > theFilledMins;
  std::map< G4int, std::map< G4int, G4int > > theFilledMaxs;
  G4int theNVoxels;
  G4double dimX, dimY, dimZ;
  G4double offsetX, offsetY, offsetZ;
  size_t* mateIDs;
  std::vector<G4Material*> thePhantomMaterials; 

};

#endif
