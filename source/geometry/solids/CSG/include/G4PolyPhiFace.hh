//
// G4PolyPhiFace.hh
//
// Definition of a face that bounds a polycone or polyhedra when it has a phi
// opening.
//
// Specifically: a face that lies on a plane that passes through
// the z axis. It has boundaries that are straight lines of arbitrary length
// and direction, but with corners aways on the same side of the z axis.
//
#ifndef G4PolyPhiFace_hh
#define G4PolyPhiFace_hh

#include "G4VCSGface.hh"

typedef struct {
        G4double r, z;          // position
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

class G4PolyPhiFace : public G4VCSGface {

	public:
	G4PolyPhiFace( const G4double *r, const G4double *z, const G4int n,
		       const G4double phi, const G4double deltaPhi, const G4bool start );
	virtual ~G4PolyPhiFace();
	
	G4bool Intersect( const G4ThreeVector &p, const G4ThreeVector &v,
			  const G4bool outgoing, const G4double surfTolerance,
			  G4double &distance, G4double &distFromSurface,
			  G4ThreeVector &normal, G4bool &allBehind );

        G4double Distance( const G4ThreeVector &p, const G4bool outgoing );
	
	EInside Inside( const G4ThreeVector &p, const G4double tolerance, 
			G4double *bestDistance );
		
	G4ThreeVector Normal( const G4ThreeVector &p,  G4double *bestDistance );

	G4double Extent( const G4ThreeVector axis );
	
	void CalculateExtent( const EAxis axis, 
			      const G4VoxelLimits &voxelLimit,
			      const G4AffineTransform &tranform,
			      G4double &min, G4double &max        );
	
	protected:
	G4PolyPhiFaceEdge	*edges;		// The edges of the face
        G4PolyPhiFaceVertex     *corners;       // And the corners
	G4int			numEdges;	// Number of edges
	G4ThreeVector		normal;		// Normal unit vector
	G4ThreeVector		radial;		// Unit vector along radial direction
	G4ThreeVector		surface;	// Point on surface
	G4double		rMin, rMax,	// Extent in r
				zMin, zMax;	// Extent in z
				
	G4bool InsideEdges( const G4double r, const G4double z );
	G4bool InsideEdges( const G4double r, const G4double z, 
			    G4double *distRZ2, G4PolyPhiFaceVertex **base3Dnorm=0,
			    G4ThreeVector **head3Dnorm=0 );
};

#endif
