//
// G4ClippablePolygon.hh
//
// Declaration of a utility class of a polycon that can be clipped
// by a voxel
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
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

	void PartialClip( const G4VoxelLimits &voxelLimit, const EAxis IgnoreMe );
	
	void ClipAlongOneAxis( const G4VoxelLimits &voxelLimit, const EAxis axis );
	
	G4bool GetExtent( const EAxis axis, 
			  G4double &min, G4double &max );
			       
	G4bool GetNumVertices() const { return vertices.entries(); }
			      
	private:
	G4ThreeVectorList vertices;
	
	void ClipToSimpleLimits( G4ThreeVectorList& pPolygon,
				 G4ThreeVectorList& outputPolygon,
				 const G4VoxelLimits& pVoxelLimit  );

};

#endif
