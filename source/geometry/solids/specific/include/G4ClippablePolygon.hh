// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ClippablePolygon.hh,v 1.4 2000-11-02 16:54:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// --------------------------------------------------------------------
// GEANT 4 class header file
//
// G4ClippablePolygon
//
// Class description:
//
//   Declaration of a utility class of a polygon that can be
//   clipped by a voxel.

// Author:
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4ClippablePolygon_hh
#define G4ClippablePolygon_hh

#include "globals.hh"
#include "geomdefs.hh"
#include "G4ThreeVector.hh"

class G4VoxelLimits;
class G4AffineTransform;

class G4AffineTransform;
class G4VoxelLimits;

#include "g4rw/tvordvec.h"
typedef G4RWTValOrderedVector<G4ThreeVector> G4ThreeVectorList;

class G4ClippablePolygon
{
  public:  // with description

	G4ClippablePolygon();
	virtual ~G4ClippablePolygon();
          // Constructor & virtual destructor.
	
	virtual void AddVertexInOrder( const G4ThreeVector vertex );
	virtual void ClearAllVertices();
	
	inline void SetNormal( const G4ThreeVector &newNormal );
	inline const G4ThreeVector GetNormal() const;
	
	virtual G4bool Clip( const G4VoxelLimits &voxelLimit );

	virtual G4bool PartialClip( const G4VoxelLimits &voxelLimit,
	                            const EAxis IgnoreMe );
          // Clip, while ignoring the indicated axis.

	virtual void ClipAlongOneAxis( const G4VoxelLimits &voxelLimit,
	                               const EAxis axis );
          // Clip along just one axis, as specified in voxelLimit.

	virtual G4bool GetExtent( const EAxis axis, 
		  		  G4double &min, G4double &max ) const;

	virtual const G4ThreeVector *GetMinPoint( const EAxis axis ) const;
          // Returns pointer to minimum point along the specified axis.
          // Take care! Do not use pointer after destroying parent polygon.

	virtual const G4ThreeVector *GetMaxPoint( const EAxis axis ) const;
          // Returns pointer to maximum point along the specified axis.
          // Take care! Do not use pointer after destroying parent polygon.

	inline G4int GetNumVertices() const;
	inline G4bool Empty() const;
	
	virtual G4bool InFrontOf( const G4ClippablePolygon &other, EAxis axis ) const;
          // Decide if the polygon is in "front" of another when
          // viewed along the specified axis. For our purposes here,
          // it is sufficient to use the minimum extent of the
          // polygon along the axis to determine this.

	virtual G4bool BehindOf( const G4ClippablePolygon &other, EAxis axis ) const;
          // Decide if this polygon is behind another.
          // Remarks in method "InFrontOf" are valid here too.

	virtual G4bool GetPlanerExtent( const G4ThreeVector &pointOnPlane, 
		  			const G4ThreeVector &planeNormal,
		  			G4double &min, G4double &max ) const;
          // Get min/max distance in or out of a plane.

  protected:  // with description

	void ClipToSimpleLimits( G4ThreeVectorList& pPolygon,
				 G4ThreeVectorList& outputPolygon,
				 const G4VoxelLimits& pVoxelLimit  );
          // pVoxelLimits must be only limited along one axis, and either
          // the maximum along the axis must be +kInfinity, or the minimum
          // -kInfinity

  protected:

	G4ThreeVectorList vertices;
	G4ThreeVector normal;
};

#include "G4ClippablePolygon.icc"

#endif
