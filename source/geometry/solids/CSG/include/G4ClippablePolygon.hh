//
// G4ClippablePolygon.hh
//
// Declaration of a utility class of a polycon that can be clipped
// by a voxel
//
#ifndef G4ClippablePolycon_hh
#define G4ClippablePolycon_hh

#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4VoxelLimits;
class G4AffineTransform;

class G4AffineTransform;
class G4VoxelLimits;

#include <rw/tvordvec.h>
typedef RWTValOrderedVector<G4ThreeVector> G4ThreeVectorList;

class G4ClippablePolygon {
	public:
	G4ClippablePolygon() {;}
	~G4ClippablePolygon() {;}
	
	void AddVertexInOrder( const G4ThreeVector vertex );
	void ClearAllVertices();
	
	void Clip( const G4VoxelLimits &voxelLimit );
	
	void GetExtent( const EAxis axis, 
			G4double &min, G4double &max );
			      
	private:
	G4ThreeVectorList vertices;
	
	void ClipToSimpleLimits( G4ThreeVectorList& pPolygon,
				 G4ThreeVectorList& outputPolygon,
				 const G4VoxelLimits& pVoxelLimit  );

};

#endif
