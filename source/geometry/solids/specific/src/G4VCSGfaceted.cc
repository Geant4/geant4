// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VCSGfaceted.cc,v 1.2 2000-04-11 16:03:40 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4VCSGfaceted.cc
//
// Implementation of the virtual class of a CSG type shape that is built
// entirely out of G4VCSGface faces.
//
// --------------------------------------------------------------------

#include "G4VCSGfaceted.hh"
#include "G4VCSGface.hh"
#include "G4SolidExtentList.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4Polyhedron.hh"   
#include "G4VGraphicsScene.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"

//
// Destructor
//
G4VCSGfaceted::~G4VCSGfaceted()
{
	DeleteStuff();
}


//
// Copy constructor
//
G4VCSGfaceted::G4VCSGfaceted( const G4VCSGfaceted &source ) : G4CSGSolid( source )
{
	CopyStuff( source );
}


//
// Assignment operator
//
const G4VCSGfaceted &G4VCSGfaceted::operator=( const G4VCSGfaceted &source )
{
	if (&source == this) return *this;
	
	DeleteStuff();
	CopyStuff( source );
	
	return *this;
}


//
// CopyStuff (protected)
//
// Copy the contents of source
//
void G4VCSGfaceted::CopyStuff( const G4VCSGfaceted &source )
{
	numFace = source.numFace;
	if (numFace == 0) return;		// odd, but permissable?
	
	faces = new G4VCSGface*[numFace];
	
	G4VCSGface **face = faces,
		   **sourceFace = source.faces;
	do {
		*face = (*sourceFace)->Clone();
	} while( ++sourceFace, ++face < faces+numFace );
}


//
// DeleteStuff (protected)
//
// Delete all allocated objects
//
void G4VCSGfaceted::DeleteStuff()
{
	if (numFace) {
		G4VCSGface **face = faces;
		do {
			delete *face;
		} while( ++face < faces + numFace );

		delete [] faces;
	}
}


//
// CalculateExtent
//
G4bool G4VCSGfaceted::CalculateExtent( const EAxis axis,
				       const G4VoxelLimits &voxelLimit,
				       const G4AffineTransform &transform,
				       G4double &min, G4double &max ) const
{
	G4SolidExtentList	extentList( axis, voxelLimit );

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
	return extentList.GetExtent( min, max );
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
	G4VCSGface *bestFace;
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
				bestFace = *face;
				if (distFromSurface <= 0) return 0;
			}
		}
	} while( ++face < faces + numFace );
	
	if (distance < kInfinity && distFromSurface<kCarTolerance/2) {
		if (bestFace->Distance(p,false) < kCarTolerance/2) distance = 0;
	}

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
	G4VCSGface *bestFace;
	
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
				bestFace = *face;
				if (distFromSurface <= 0) break;
			}
		}
	} while( ++face < faces + numFace );
	
	if (distance < kInfinity) {
		if (distFromSurface <= 0)
			distance = 0;
		else if (distFromSurface<kCarTolerance/2) {
			if (bestFace->Distance(p,true) < kCarTolerance/2) distance = 0;
		}

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
