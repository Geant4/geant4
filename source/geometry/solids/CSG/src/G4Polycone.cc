//
// G4Polycone.cc
//
// Implementation of a CSG polycone
//
#include "G4Polycone.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyPhiFace.hh"

#include "G4Polyhedron.hh"


//
// Constructor (GEANT3 style parameters)
//	
G4Polycone::G4Polycone( G4String name, 
                        const G4double phiStart,
                        const G4double phiTotal,
                        const G4int numZPlanes,
                        const G4double zPlane[],
                        const G4double rInner[],
                        const G4double rOuter[]  ) : G4VCSGfaceted( name )
{
	//
	// Real ugly
	//
	original_parameters.exist = true;
	
	original_parameters.Start_angle = phiStart;
	original_parameters.Opening_angle = phiTotal;
	original_parameters.Num_z_planes = numZPlanes;
	original_parameters.Z_values = new G4double[numZPlanes];
	original_parameters.Rmin = new G4double[numZPlanes];
	original_parameters.Rmax = new G4double[numZPlanes];
        G4int i;
	for (i=0; i<numZPlanes; i++) {
		original_parameters.Z_values[i] = zPlane[i];
		original_parameters.Rmin[i] = rInner[i];
		original_parameters.Rmax[i] = rOuter[i];
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
		*rOut = rOuter[i];
		*rIn  = rInner[i];
		*zOut = *zIn = zPlane[i];
	}
	
	Create( phiStart, phiTotal, numZPlanes*2, r, z );
	
	delete [] r;
	delete [] z;
}


//
// Constructor (generic parameters)
//
G4Polycone::G4Polycone( G4String name, 
		    	const G4double phiStart,
                    	const G4double phiTotal,
		    	const G4int    numRZ,
		    	const G4double r[],
		    	const G4double z[]	 ) : G4VCSGfaceted( name )
{
	original_parameters.exist = false;
	
	Create( phiStart, phiTotal, numRZ, r, z );
}


//
// Create
//
// Generic create routine, called by each constructor after conversion of arguments
//
void G4Polycone::Create( const G4double phiStart,
            	     	 const G4double phiTotal,
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
	// Allocate corner array. We may not end up using all of this array,
	// since we delete duplicate corners, but that's not so bad
	//
	corners = new G4PolyconeSideRZ[numRZ];

	//
	// Copy corners, avoiding duplicates on the way
	//
	// We should also look for divided conical surfaces...
	// We must also look for overlapping surfaces...
	//
	G4PolyconeSideRZ *next = corners;
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
	// Construct conical faces
	//
	// But! Don't construct a face if both points are at zero radius!
	//
	G4PolyconeSideRZ *corner = corners,
			 *prev = corners + numCorner-1,
			 *nextNext;
	G4VCSGface	**face = faces;
	do {
		next = corner+1;
		if (next >= corners+numCorner) next = corners;
		nextNext = next+1;
		if (nextNext >= corners+numCorner) nextNext = corners;
		
		if (corner->r < 1/kInfinity && next->r < 1/kInfinity) continue;
		
		*face++ = new G4PolyconeSide( prev, corner, next, nextNext,
					      startPhi, endPhi-startPhi, phiIsOpen );
	} while( prev=corner, corner=next, corner > corners );
	
	if (phiIsOpen) {
		//
		// Construct phi open edges
		//
		*face++ = new G4PolyPhiFace( r, z, numRZ, startPhi, 0, true  );
		*face++ = new G4PolyPhiFace( r, z, numRZ, endPhi,   0, false );
	}
	
	//
	// We might have dropped a face or two: recalculate numFace
	//
	numFace = face-faces;
}


//
// Destructor
//
G4Polycone::~G4Polycone()
{
	delete [] corners;
	
	if (original_parameters.exist) {
		delete [] original_parameters.Z_values;
		delete [] original_parameters.Rmin;
		delete [] original_parameters.Rmax;
	}
}


//
// ComputeDimensions
//
void G4Polycone::ComputeDimensions( G4VPVParameterisation* p,
				    const G4int n,
				    const G4VPhysicalVolume* pRep)
{
}


//
// CreatePolyhedron
//
G4Polyhedron *G4Polycone::CreatePolyhedron() const
{ 
	//
	// It is *really* unfortunate how the design in /graphics_reps is
	// written to parallel the design in /geometry/solids. Ugly, ugly, ugly.
	//
	// This has to be fixed, but I won't do it now. Fake it for the moment.
	// 
	if (original_parameters.exist) {
	
		return new G4PolyhedronPcon( original_parameters.Start_angle,
					     original_parameters.Opening_angle,
					     original_parameters.Num_z_planes,
					     original_parameters.Z_values,
					     original_parameters.Rmin,
					     original_parameters.Rmax);
	}
	else {
		G4Exception( "G4Polycone: waiting for graphics_reps to catch up" );
		return 0;
	}
}	


//
// CreateNURBS
//
G4NURBS *G4Polycone::CreateNURBS() const
{
	return 0;
}
