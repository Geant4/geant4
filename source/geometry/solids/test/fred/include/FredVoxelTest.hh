//
// FredVoxelTest.hh
//
// Declaration of Fred's voxel solid tester
//

#ifndef FredVoxelTest_hh
#define FredVoxelTest_hh

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

class G4VSolid;
class G4VVisManager;

class FredVoxelTest {
	public:
	FredVoxelTest(  );
	~FredVoxelTest();
	
	void Test( const EAxis axis, const G4VSolid *solid );
	void Draw();
	
	void SetExtent( const EAxis axis, const G4double min, const G4double max );
	void SetOrigin( const G4ThreeVector origin );
	void Rotate( const EAxis axis, const G4double value );
	void ResetRotation();
	
	protected:
	
	void PlotMarker( G4ThreeVector point, 
			 G4VVisManager *visManager, G4bool test=false );
	void PlotLine( G4ThreeVector start, G4ThreeVector end, 
		       G4VVisManager *visManager,G4bool test=false );
	
	G4AffineTransform transform,
			  inverseTransform;
	G4VoxelLimits	  voxelLimits;
	G4RotationMatrix  rotation;
	
	typedef struct {
		const G4VSolid *solid;
		EAxis	 axis;
		G4bool	 result;
		G4double min, max;
	} FredVoxelTestResult;
	
	FredVoxelTestResult test;
};	

#endif
