//
// G4VCSGfaceted.cc
//
// Implementation of the virtual class of a CSG type shape that is built
// entirely out of G4VCSGface faces.
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#include "G4VCSGfaceted.hh"
#include "G4VCSGface.hh"
#include "G4SolidExtentList.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4Polyhedron.hh"   
#include "G4VGraphicsScene.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

//
// Destructor
//
G4VCSGfaceted::~G4VCSGfaceted()
{
	G4VCSGface **face = faces;
	do {
		delete *face;
	} while( ++face < faces + numFace );
	
	delete [] faces;
}


//
// CalculateExtent
//
G4bool G4VCSGfaceted::CalculateExtent( const EAxis axis,
				       const G4VoxelLimits &voxelLimit,
				       const G4AffineTransform &transform,
				       G4double &min, G4double &max ) const
{
	G4SolidExtentList	extentList;

	//
	// Loop over all faces, checking min/max extent as we go.
	//
	G4VCSGface **face = faces;
	do {
		(*face)->CalculateExtent( axis, voxelLimit, transform, extentList );
	} while( ++face < faces + numFace );
	
	//
	// Return min/max value
	//
	G4bool answer;
	
	if (voxelLimit.IsLimited(axis)) {
		max = voxelLimit.GetMaxExtent(axis);
		min = voxelLimit.GetMinExtent(axis);
		answer = extentList.UpdateLimitedExtent( min, max );
	}
	else 
		answer = extentList.GetExtent( min, max );
		
	return answer;
}


//
// Inside
//
// It could be a good idea to override this virtual
// member to add first a simple test (such as spherical
// test or whatnot) and to call this version only if
// the simplier test fails.
//
EInside G4VCSGfaceted::Inside( const G4ThreeVector &p ) const
{
	EInside answer;
	G4VCSGface **face = faces;
	G4double best = kInfinity;
	do {
		G4double distance;
		EInside result = (*face)->Inside( p, kCarTolerance/2, &distance );
		if (result == kSurface) return kSurface;
		if (distance < best) {
			best = distance;
			answer = result;
		}
	} while( ++face < faces + numFace );

	return answer;
}


//
// SurfaceNormal
//
G4ThreeVector G4VCSGfaceted::SurfaceNormal( const G4ThreeVector& p) const
{
	G4ThreeVector answer;
	G4VCSGface **face = faces;
	G4double best = kInfinity;
	do {
		G4double distance;
		G4ThreeVector normal = (*face)->Normal( p, &distance );
		if (distance < best) {
			best = distance;
			answer = normal;
		}
	} while( ++face < faces + numFace );

	return answer;
}


//
// DistanceToIn(p,v)
//
G4double G4VCSGfaceted::DistanceToIn( const G4ThreeVector &p, const G4ThreeVector &v ) const
{
	G4double distance = kInfinity;
	G4double distFromSurface;
	G4VCSGface **face = faces;
	do {
		G4double 	faceDistance,
			 	faceDistFromSurface;
		G4ThreeVector 	faceNormal;
		G4bool		faceAllBehind;
		if ((*face)->Intersect( p, v, false, kCarTolerance/2,
				        faceDistance, faceDistFromSurface,
				        faceNormal, faceAllBehind ) ) {
			//
			// Intersecting face
			//
			if (faceDistance < distance) {
				distance = faceDistance;
				distFromSurface = faceDistFromSurface;
			}
		}
	} while( ++face < faces + numFace );
	
	if ((distance < kInfinity) && (fabs(distFromSurface)<kCarTolerance/2) ) distance = 0;

	return distance;
}


//
// DistanceToIn(p)
//
G4double G4VCSGfaceted::DistanceToIn( const G4ThreeVector &p ) const
{
	return DistanceTo( p, false );
}


//
// DistanceToOut(p,v)
//
G4double G4VCSGfaceted::DistanceToOut( const G4ThreeVector &p, const G4ThreeVector &v,
				       const G4bool calcNorm,
				       G4bool *validNorm, G4ThreeVector *n ) const
{
	G4bool allBehind = true;
	G4double distance = kInfinity;
	G4double distFromSurface;
	G4ThreeVector normal;
	
	G4VCSGface **face = faces;
	do {
		G4double	faceDistance,
				faceDistFromSurface;
		G4ThreeVector	faceNormal;
		G4bool		faceAllBehind;
		if ((*face)->Intersect( p, v, true, kCarTolerance/2,
				        faceDistance, faceDistFromSurface,
				        faceNormal, faceAllBehind ) ) {
			//
			// Intersecting face
			//
			if ( (distance < kInfinity) || (!faceAllBehind) ) allBehind = false;
			if (faceDistance < distance) {
				distance = faceDistance;
				distFromSurface = faceDistFromSurface;
				normal = faceNormal;
			}
		}
	} while( ++face < faces + numFace );
	
	if (distance < kInfinity) {
		if (fabs(distFromSurface)<kCarTolerance/2) distance = 0;

		if (calcNorm) {
			*validNorm = allBehind;
			*n = normal;
		}
	}
	else {
		if (calcNorm) *validNorm = false;
	}

	return distance;
}


//
// DistanceToOut(p)
//
G4double G4VCSGfaceted::DistanceToOut( const G4ThreeVector &p ) const
{
	return DistanceTo( p, true );
}


//
// DistanceTo
//
// Protected routine called by DistanceToIn and DistanceToOut
//
G4double G4VCSGfaceted::DistanceTo( const G4ThreeVector &p, const G4bool outgoing ) const
{
	G4VCSGface **face = faces;
	G4double best = kInfinity;
	do {
		G4double distance = (*face)->Distance( p, outgoing );
		if (distance < best) best = distance;
	} while( ++face < faces + numFace );

	return (best < 0.5*kCarTolerance) ? 0 : best;
}


//
// DescribeYourselfTo
//
void G4VCSGfaceted::DescribeYourselfTo( G4VGraphicsScene& scene ) const
{
   scene.AddThis( *this );
}


//
// GetExtent
//
G4VisExtent G4VCSGfaceted::GetExtent() const
{  
	G4ThreeVector plusX(1,0,0), minusX(-1,0,0),
		      plusY(0,1,0), minusY(0,-1,0),
		      plusZ(0,0,1), minusZ(0,0,-0);
	G4double answer, maxX = -kInfinity, 
			 minX = -kInfinity, 
			 maxY = -kInfinity, 
			 minY = -kInfinity,
			 maxZ = -kInfinity, 
			 minZ = -kInfinity;
	
	//
	// Ask everyone about x, y, and z
	//
	G4VCSGface **face = faces;
	do {
		answer = (*face)->Extent( plusX );
		if (answer > maxX) maxX = answer;

		answer = (*face)->Extent( minusX );
		if (answer > minX) minX = answer;

		answer = (*face)->Extent( plusY );
		if (answer > maxY) maxY = answer;

		answer = (*face)->Extent( minusY );
		if (answer > minY) minY = answer;

		answer = (*face)->Extent( plusZ );
		if (answer > maxZ) maxZ = answer;

		answer = (*face)->Extent( minusZ );
		if (answer > minZ) minZ = answer;
	} while( ++face < faces + numFace );
	
	//
	// Yes, the signs are correct!
	//
	return G4VisExtent( -minX, maxY, -minY, maxY, -minZ, maxZ );
}
