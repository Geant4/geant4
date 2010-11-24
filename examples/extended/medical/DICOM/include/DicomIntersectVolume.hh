#ifndef DicomIntersectVolume__HH
#define DicomIntersectVolume__HH

#include "G4UImessenger.hh"

#include <vector>
#include <fstream>
#include "G4PhantomParameterisation.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PhantomParameterisation;
class G4UIcmdWithAString;
class G4VSolid;

class DicomIntersectVolume : public G4UImessenger
{
public:
  DicomIntersectVolume();
  ~DicomIntersectVolume();

  void SetNewValue(G4UIcommand * command,
		   G4String newValues);

private:

  void BuildUserSolid( std::vector<G4String> params );
  void BuildG4Solid( std::vector<G4String> params );
  G4PhantomParameterisation* GetPhantomParam(G4bool bMustExist);
  G4bool IsPhantomVolume( G4VPhysicalVolume* pv );
  std::vector<G4VPhysicalVolume*> GetPhysicalVolumes( const G4String& name, bool exists, G4int nVols );
  std::vector<G4LogicalVolume*> GetLogicalVolumes( const G4String& name, bool exists, G4int nVols );
  std::vector<G4String> GetWordsInString( const G4String& stemp);

private:
  G4UIcmdWithAString* theUserVolumeCmd;
  G4UIcmdWithAString* theG4VolumeCmd;

  G4VSolid* theSolid;

  G4VPhysicalVolume* thePhantomVolume;

  std::ofstream fout;

  G4PhantomParameterisation* theRegularParam;

  G4ThreeVector thePhantomMinusCorner;

  G4bool* theVoxelIsInside;
};

#endif
