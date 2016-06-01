//
// G4Polyhedra.cc
//
// Implementation of a CSG polyhedra, as an inherited class of G4VCSGfaceted.
//
// To be done:
//    * Checks for bad input should be improved. It is now possible for
//      users to specify crazy polyhedra parameters without complaint that
//      could produce unpredictable results.
//    * Cracks: there are probably small cracks in the seams between the
//      phi face (G4PolyPhiFace) and sides (G4PolyhedraSide) that are not
//      entirely leakproof. Also, I am not sure all vertices are leak proof.
//    * Many optimizations are possible, but not implemented.
//    * Visualization needs to be updated outside of this routine.
//
#include "G4Polyhedra.hh"
#include "G4PolyhedraSide.hh"
#include "G4PolyPhiFace.hh"

#include "G4Polyhedron.hh"
#include "G4EnclosingCylinder.hh"


//
// Constructor (GEANT3 style parameters)
//
// GEANT3 PGON radii are specified in the distance to the norm of each face.
//	
G4Polyhedra::G4Polyhedra( G4String name, 
                          const G4double phiStart,
                          const G4double thePhiTotal,
 		          const G4double theNumSide,	
                          const G4int numZPlanes,
                          const G4double zPlane[],
                          const G4double rInner[],
                          const G4double rOuter[]  ) : G4VCSGfaceted( name )
{
	if (theNumSide <= 0) G4Exception( "G4Polyhedra:: must have at least one side" );

	//
	// Calculate conversion factor from G3 radius to G4 radius
	//
	G4double phiTotal = thePhiTotal;
	if (phiTotal <=0 || phiTotal > 2*M_PI) phiTotal = 2*M_PI;
	G4double convertRad = cos(0.5*phiTotal/theNumSide);

	//
	// Real ugly
	//
	original_parameters.exist = true;
	
	original_parameters.numSide = theNumSide;
	original_parameters.Start_angle = phiStart;
	original_parameters.Opening_angle = phiTotal;
	original_parameters.Num_z_planes = numZPlanes;
	original_parameters.Z_values = new G4double[numZPlanes];
	original_parameters.Rmin = new G4double[numZPlanes];
	original_parameters.Rmax = new G4double[numZPlanes];

        G4int i;
	for (i=0; i<numZPlanes; i++) {
		original_parameters.Z_values[i] = zPlane[i];
		original_parameters.Rmin[i] = rInner[i]/convertRad;
		original_parameters.Rmax[i] = rOuter[i]/convertRad;
	}
		
	//
	// Translate GEANT3 into generic parameters
	// Duplicate vertices and divided surfaces are (or should be) dealt with
	// by routine "Create."
	//
	G4double *r = new G4double[numZPlanes*2];
	G4double *z = new G4double[numZPlanes*2];
	
	G4double *rOut = r + numZPlanes,
		 *zOut = z + numZPlanes,
		 *rIn = rOut-1,
		 *zIn = zOut-1;
		 
	for( i=0; i < numZPlanes; i++, rOut++, zOut++, rIn--, zIn-- ) {
		*rOut = rOuter[i]/convertRad;
		*rIn  = rInner[i]/convertRad;
		*zOut = *zIn = zPlane[i];
	}
	
	Create( phiStart, phiTotal, theNumSide, numZPlanes*2, r, z );
	
	delete [] r;
	delete [] z;
}


//
// Constructor (generic parameters)
//
G4Polyhedra::G4Polyhedra( G4String name, 
		    	  const G4double phiStart,
                    	  const G4double phiTotal,
 		          const G4double theNumSide,	
		    	  const G4int    numRZ,
		    	  const G4double r[],
		    	  const G4double z[]	 ) : G4VCSGfaceted( name )
{
	original_parameters.exist = false;
	
	Create( phiStart, phiTotal, theNumSide, numRZ, r, z );
}


//
// Create
//
// Generic create routine, called by each constructor after conversion of arguments
//
void G4Polyhedra::Create( const G4double phiStart,
            	     	  const G4double phiTotal,
 		          const G4double theNumSide,	
		     	  const G4int	numRZ,
		     	  const G4double r[],
	             	  const G4double z[]	  )
{
        //
        // Phi opening? Account for some possible roundoff, and interpret
        // nonsense value as representing no phi opening
        //
        if (phiTotal <= 0 || phiTotal > 2.0*M_PI-1E-10) {
                phiIsOpen = false;
		startPhi = 0;
                endPhi = 2*M_PI;
        }
        else {
                phiIsOpen = true;
		
		//
		// Convert phi into our convention
		//
		startPhi = phiStart;
		while( startPhi < 0 ) startPhi += 2*M_PI;
		
		endPhi = phiStart+phiTotal;
		while( endPhi < startPhi ) endPhi += 2*M_PI;
        }
	
	//
	// Save number sides
	//
	numSide = theNumSide;
	
	//
	// Allocate corner array. We may not end up using all of this array,
	// since we delete duplicate corners, but that's not so bad
	//
	corners = new G4PolyhedraSideRZ[numRZ];

	//
	// Copy corners, avoiding duplicates on the way
	//
	// We should also look for divided conical surfaces...
	// We must also look for overlapping surfaces...
	//
	G4PolyhedraSideRZ *next = corners;
	const G4double *rOne = r;
	const G4double *zOne = z;
	const G4double *rNext, *zNext;
	G4bool	 notFinished;
	do {
		rNext = rOne + 1;
		zNext = zOne + 1;
		if (notFinished = (rNext < r+numRZ)) {
			if (*rNext == *rOne && *zNext == *zOne) continue;
		}
		
		next->r = *rOne;
		next->z = *zOne;
		next++;
	} while( rOne=rNext, zOne=zNext, notFinished );
			
	numCorner = next - corners;
	
	//
	// Allocate face pointer array
	//
	numFace = phiIsOpen ? numCorner+2 : numCorner;
	faces = new G4VCSGface*[numFace];
	
	//
	// Construct side faces
	//
	// To do so properly, we need to keep track of four successive RZ
	// corners.
	//
	// But! Don't construct a face if both points are at zero radius!
	//
	G4PolyhedraSideRZ *corner = corners,
			  *prev = corners + numCorner-1,
			  *nextNext;
	G4VCSGface	 **face = faces;
	do {
		next = corner+1;
		if (next >= corners+numCorner) next = corners;
		nextNext = next+1;
		if (nextNext >= corners+numCorner) nextNext = corners;
		
		if (corner->r < 1/kInfinity && next->r < 1/kInfinity) continue;
		
		*face++ = new G4PolyhedraSide( prev, corner, next, nextNext,
					       numSide, startPhi, endPhi-startPhi, phiIsOpen );
	} while( prev=corner, corner=next, corner > corners );
	
	if (phiIsOpen) {
		//
		// Construct phi open edges
		//
		*face++ = new G4PolyPhiFace( r, z, numRZ, startPhi, phiTotal/numSide, true  );
		*face++ = new G4PolyPhiFace( r, z, numRZ, endPhi,   phiTotal/numSide, false );
	}
	
	//
	// We might have dropped a face or two: recalculate numFace
	//
	numFace = face-faces;
	
	//
	// Make enclosingCylinder
	//
	enclosingCylinder = new G4EnclosingCylinder( r, z, numRZ, 
						     phiIsOpen, phiStart, phiTotal );
}


//
// Destructor
//
G4Polyhedra::~G4Polyhedra()
{
	delete [] corners;

	if (original_parameters.exist) {
		delete [] original_parameters.Z_values;
		delete [] original_parameters.Rmin;
		delete [] original_parameters.Rmax;
	}
	
	delete enclosingCylinder;
}


//
// Inside
//
EInside G4Polyhedra::Inside( const G4ThreeVector &p ) const
{
	//
	// Quick test
	//
	if (enclosingCylinder->Outside(p)) return kOutside;

	//
	// Long answer
	//
	return G4VCSGfaceted::Inside(p);
}


//
// DistanceToIn
//
G4double G4Polyhedra::DistanceToIn( const G4ThreeVector &p, const G4ThreeVector &v ) const
{
	//
	// Quick test
	//
	if (enclosingCylinder->Misses(p,v)) return kInfinity;
	
	//
	// Long answer
	//
	return G4VCSGfaceted::DistanceToIn( p, v );
}


//
// ComputeDimensions
//
void G4Polyhedra::ComputeDimensions( G4VPVParameterisation* p,
				    const G4int n,
				    const G4VPhysicalVolume* pRep)
{
}


//
// CreatePolyhedron
//
G4Polyhedron *G4Polyhedra::CreatePolyhedron() const
{ 
	//
	// It is *really* unfortunate how the design in /graphics_reps is
	// written to parallel the design in /geometry/solids. Ugly, ugly, ugly.
	//
	// This has to be fixed, but I won't do it now. Fake it for the moment.
	// 
	if (original_parameters.exist) {
	
		return new G4PolyhedronPgon( original_parameters.Start_angle,
					     original_parameters.Opening_angle,
					     original_parameters.numSide,
					     original_parameters.Num_z_planes,
					     original_parameters.Z_values,
					     original_parameters.Rmin,
					     original_parameters.Rmax);
	}
	else {
		G4Exception( "G4Polyhedra: waiting for graphics_reps to catch up" );
		return 0;
	}

}	


//
// CreateNURBS
//
G4NURBS *G4Polyhedra::CreateNURBS() const
{
	return 0;
}
