// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PolyPhiFace.hh,v 1.3 2000-11-02 16:54:48 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class header file
//
//
// G4PolyPhiFace
//
// Class description:
//
//   Definition of a face that bounds a polycone or polyhedra when
//   it has a phi opening:
//
//   G4PolyPhiFace( const G4ReduciblePolygon *rz,
//                        G4double phi,
//                        G4double deltaPhi,
//                        G4double phiOther )
//
//   Specifically: a face that lies on a plane that passes through
//   the z axis. It has boundaries that are straight lines of arbitrary
//   length and direction, but with corners aways on the same side of
//   the z axis.

// Author: 
//   David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#ifndef G4PolyPhiFace_hh
#define G4PolyPhiFace_hh

#include "G4VCSGface.hh"

class G4ReduciblePolygon;

typedef struct {
        G4double x, y, r, z;          // position
        G4double rNorm, 
                 zNorm;         // r/z normal
	G4ThreeVector norm3D;	// 3D normal
} G4PolyPhiFaceVertex;

typedef struct {
        G4PolyPhiFaceVertex     *v0, *v1;       // Corners
        G4double tr, tz,                        // Unit vector along edge
                 length;                        // Length of edge
	G4ThreeVector norm3D;			// 3D edge normal vector
} G4PolyPhiFaceEdge;

class G4PolyPhiFace : public G4VCSGface
{

  public:  // with description

        G4PolyPhiFace( const G4ReduciblePolygon *rz,
                       G4double phi, G4double deltaPhi, G4double phiOther );
          // Constructor.
          // Points r,z should be supplied in clockwise order in r,z.
          // For example:
          //                [1]---------[2]         ^ R
          //                 |           |          |
          //                 |           |          +--> z
          //                [0]---------[3]

        virtual ~G4PolyPhiFace();
          // Destructor. Removes edges and corners.

	G4PolyPhiFace( const G4PolyPhiFace &source );
	G4PolyPhiFace& operator=( const G4PolyPhiFace &source );
	  // Copy constructor and assgnment operator.

	G4bool Intersect( const G4ThreeVector &p, const G4ThreeVector &v,
			  G4bool outgoing, G4double surfTolerance,
			  G4double &distance, G4double &distFromSurface,
			  G4ThreeVector &normal, G4bool &allBehind );

        G4double Distance( const G4ThreeVector &p, G4bool outgoing );
	
	EInside Inside( const G4ThreeVector &p, G4double tolerance, 
			G4double *bestDistance );
		
	G4ThreeVector Normal( const G4ThreeVector &p,  G4double *bestDistance );

	G4double Extent( const G4ThreeVector axis );
	
	void CalculateExtent( const EAxis axis, 
			      const G4VoxelLimits &voxelLimit,
			      const G4AffineTransform &tranform,
			      G4SolidExtentList &extentList        );

	inline G4VCSGface *Clone();
          // Allocates on the heap a clone of this face.

  public:  // without description

	void Diagnose( G4VSolid *solid );
          // Throw an exception if something is found inconsistent with
          // the solid. For debugging purposes only


  protected:

	G4PolyPhiFaceEdge	*edges;	    // The edges of the face
        G4PolyPhiFaceVertex     *corners;   // And the corners
	G4int			numEdges;   // Number of edges
	G4ThreeVector		normal;     // Normal unit vector
	G4ThreeVector		radial;	    // Unit vector along radial direction
	G4ThreeVector		surface;    // Point on surface
	G4double		rMin, rMax, // Extent in r
				zMin, zMax; // Extent in z
	G4bool			allBehind;  // True if the polycone/polyhedra
                                            // is behind the place of this face

        G4bool InsideEdgesExact( G4double r, G4double z, G4double normSign,
                                 const G4ThreeVector &p, const G4ThreeVector &v );
          // Decide if the point in r,z is inside the edges of our face,
          // **but** do so consistently with other faces.

	G4bool InsideEdges( G4double r, G4double z );
	G4bool InsideEdges( G4double r, G4double z, G4double *distRZ2,
                            G4PolyPhiFaceVertex **base3Dnorm=0,
			    G4ThreeVector **head3Dnorm=0 );
          // Decide if the point in r,z is inside the edges of our face.

	inline G4double ExactZOrder( G4double z, 
				     G4double qx, G4double qy, G4double qz, 
				     const G4ThreeVector &v, 
				     G4double normSign,
				     const G4PolyPhiFaceVertex *vert ) const;
          // Decide precisely whether a trajectory passes to the left, right,
          // or exactly passes through the z position of a vertex point in face.

	void CopyStuff( const G4PolyPhiFace &source );
};

#include "G4PolyPhiFace.icc"

#endif
