//
// G4Polyhedra.cc
//
// Implementation of a CSG polyhedra, as an inherited class of G4VCSGfaceted.
//
// To be done:
//    * Cracks: there are probably small cracks in the seams between the
//      phi face (G4PolyPhiFace) and sides (G4PolyhedraSide) that are not
//      entirely leakproof. Also, I am not sure all vertices are leak proof.
//    * Many optimizations are possible, but not implemented.
//    * Visualization needs to be updated outside of this routine.
//
// Utility classes:
//    * G4EnclosingCylinder: I decided a quick check of geometry would be a
//      good idea (for CPU speed). If the quick check fails, the regular
//      full-blown G4VCSGfaceted version is invoked.
//    * G4ReduciblePolygon: Really meant as a check of input parameters,
//      this utility class also "converts" the GEANT3-like PGON/PCON
//      arguments into the newer ones.
// Both these classes are implemented outside this file because they are
// shared with G4Polycone.
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
#include "G4Polyhedra.hh"
#include "G4PolyhedraSide.hh"
#include "G4PolyPhiFace.hh"

#include "G4Polyhedron.hh"
#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"

//
// Constructor (GEANT3 style parameters)
//
// GEANT3 PGON radii are specified in the distance to the norm of each face.
//	
G4Polyhedra::G4Polyhedra( G4String name, 
                          const G4double phiStart,
                          const G4double thePhiTotal,
 		          const G4int theNumSide,	
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
	if (phiTotal <=0 || phiTotal >= 2*M_PI*(1-DBL_EPSILON)) phiTotal = 2*M_PI;
	G4double convertRad = cos(0.5*phiTotal/theNumSide);

	//
	// Some historical stuff
	//
	original_parameters = new G4PolyhedraHistorical;
	
	original_parameters->numSide = theNumSide;
	original_parameters->Start_angle = phiStart;
	original_parameters->Opening_angle = phiTotal;
	original_parameters->Num_z_planes = numZPlanes;
	original_parameters->Z_values = new G4double[numZPlanes];
	original_parameters->Rmin = new G4double[numZPlanes];
	original_parameters->Rmax = new G4double[numZPlanes];

        G4int i;
	for (i=0; i<numZPlanes; i++) {
		original_parameters->Z_values[i] = zPlane[i];
		original_parameters->Rmin[i] = rInner[i]/convertRad;
		original_parameters->Rmax[i] = rOuter[i]/convertRad;
	}
	
	
	//
	// Build RZ polygon using special PCON/PGON GEANT3 constructor
	//
	G4ReduciblePolygon *rz = new G4ReduciblePolygon( rInner, rOuter, zPlane, numZPlanes );
	rz->ScaleA( 1/convertRad );
	
	//
	// Do the real work
	//
	Create( phiStart, phiTotal, theNumSide, rz );
	
	delete rz;
}


//
// Constructor (generic parameters)
//
G4Polyhedra::G4Polyhedra( G4String name, 
		    	  const G4double phiStart,
                    	  const G4double phiTotal,
 		          const G4int    theNumSide,	
		    	  const G4int    numRZ,
		    	  const G4double r[],
		    	  const G4double z[]	 ) : G4VCSGfaceted( name )
{
	original_parameters = 0;
	
	G4ReduciblePolygon *rz = new G4ReduciblePolygon( r, z, numRZ );
	
	Create( phiStart, phiTotal, theNumSide, rz );
	
	delete rz;
}


//
// Create
//
// Generic create routine, called by each constructor after conversion of arguments
//
void G4Polyhedra::Create( const G4double phiStart,
            	     	  const G4double phiTotal,
 		          const G4int    theNumSide,	
		     	  G4ReduciblePolygon *rz  )
{
	//
	// Perform checks of rz values
	//
	if (rz->Amin() < 0.0) 
		G4Exception( "G4Polyhedra: Illegal input parameters: All R values must be >= 0" );
		
	G4double rzArea = rz->Area();
	if (rzArea < -kCarTolerance) rz->ReverseOrder();

	else if (rzArea < -kCarTolerance)
		G4Exception( "G4Polyhedra: Illegal input parameters: R/Z cross section is zero or near zero" );
		
	if ((!rz->RemoveDuplicateVertices( kCarTolerance )) || 
	    (!rz->RemoveRedundantVertices( kCarTolerance ))     ) 
	    	G4Exception( "G4Polyhedra: Illegal input parameters: Too few unique R/Z values" );

	if (rz->CrossesItself( 1/kInfinity )) 
		G4Exception( "G4Polyhedra: Illegal input parameters: R/Z segments cross" );

	numCorner = rz->NumVertices();


	startPhi = phiStart;
	while( startPhi < 0 ) startPhi += 2*M_PI;
        //
        // Phi opening? Account for some possible roundoff, and interpret
        // nonsense value as representing no phi opening
        //
        if (phiTotal <= 0 || phiTotal > 2.0*M_PI*(1-DBL_EPSILON)) {
                phiIsOpen = false;
                endPhi = phiStart+2*M_PI;
        }
        else {
                phiIsOpen = true;
		
		//
		// Convert phi into our convention
		//
		endPhi = phiStart+phiTotal;
		while( endPhi < startPhi ) endPhi += 2*M_PI;
        }
	
	//
	// Save number sides
	//
	numSide = theNumSide;
	
	//
	// Allocate corner array. 
	//
	corners = new G4PolyhedraSideRZ[numCorner];

	//
	// Copy corners
	//
	G4ReduciblePolygonIterator iterRZ(rz);
	
	G4PolyhedraSideRZ *next = corners;
	iterRZ.Begin();
	do {
		next->r = iterRZ.GetA();
		next->z = iterRZ.GetB();
	} while( ++next, iterRZ.Next() );
	
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

		//
		// We must decide here if we can dare declare one of our faces
		// as having a "valid" normal (i.e. allBehind = true). This
		// is never possible if the face faces "inward" in r *unless*
		// we have only one side
		//
		G4bool allBehind;
		if ((corner->z > next->z) && (numSide > 1)) {
			allBehind = false;
		}
		else {
			//
			// Otherwise, it is only true if the line passing
			// through the two points of the segment do not
			// split the r/z cross section
			//
			allBehind = !rz->BisectedBy( corner->r, corner->z,
						     next->r, next->z, kCarTolerance );
		}
		
		*face++ = new G4PolyhedraSide( prev, corner, next, nextNext,
					       numSide, startPhi, endPhi-startPhi, phiIsOpen );
	} while( prev=corner, corner=next, corner > corners );
	
	if (phiIsOpen) {
		//
		// Construct phi open edges
		//
		*face++ = new G4PolyPhiFace( rz, startPhi, phiTotal/numSide, endPhi );
		*face++ = new G4PolyPhiFace( rz, endPhi,   phiTotal/numSide, startPhi );
	}
	
	//
	// We might have dropped a face or two: recalculate numFace
	//
	numFace = face-faces;
	
	//
	// Make enclosingCylinder
	//
	enclosingCylinder = new G4EnclosingCylinder( rz, phiIsOpen, phiStart, phiTotal );
}


//
// Destructor
//
G4Polyhedra::~G4Polyhedra()
{
	delete [] corners;
	if (original_parameters) delete original_parameters;
	
	delete enclosingCylinder;
}


//
// Copy constructor
//
G4Polyhedra::G4Polyhedra( const G4Polyhedra &source ) : G4VCSGfaceted( source )
{
	CopyStuff( source );
}


//
// Assignment operator
//
const G4Polyhedra &G4Polyhedra::operator=( const G4Polyhedra &source )
{
	if (this == &source) return *this;

	G4VCSGfaceted::operator=( source );
	
	delete [] corners;
	if (original_parameters) delete original_parameters;
	
	delete enclosingCylinder;
	
	CopyStuff( source );
	
	return *this;
}


//
// CopyStuff
//
void G4Polyhedra::CopyStuff( const G4Polyhedra &source )
{
	//
	// Simple stuff
	//
	numSide		= source.numSide;
	startPhi	= source.startPhi;
	endPhi		= source.endPhi;
	phiIsOpen	= source.phiIsOpen;
	numCorner	= source.numCorner;

	//
	// The corner array
	//
	corners = new G4PolyhedraSideRZ[numCorner];
	
	G4PolyhedraSideRZ	*corn = corners,
				*sourceCorn = source.corners;
	do {
		*corn = *sourceCorn;
	} while( ++sourceCorn, ++corn < corners+numCorner );
	
	//
	// Original parameters
	//
	if (source.original_parameters) {
		original_parameters = new G4PolyhedraHistorical( *source.original_parameters );
	}
	
	//
	// Enclosing cylinder
	//
	enclosingCylinder = new G4EnclosingCylinder( *source.enclosingCylinder );
}


//
// Inside
//
// This is an override of G4VCSGfaceted::Inside, created in order to speed things
// up by first checking with G4EnclosingCylinder.
//
EInside G4Polyhedra::Inside( const G4ThreeVector &p ) const
{
	//
	// Quick test
	//
	if (enclosingCylinder->MustBeOutside(p)) return kOutside;

	//
	// Long answer
	//
	return G4VCSGfaceted::Inside(p);
}


//
// DistanceToIn
//
// This is an override of G4VCSGfaceted::Inside, created in order to speed things
// up by first checking with G4EnclosingCylinder.
//
G4double G4Polyhedra::DistanceToIn( const G4ThreeVector &p, const G4ThreeVector &v ) const
{
	//
	// Quick test
	//
	if (enclosingCylinder->ShouldMiss(p,v)) return kInfinity;
	
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
	// This has to be fixed in visualization. Fake it for the moment.
	// 
	if (original_parameters) {
	
		return new G4PolyhedronPgon( original_parameters->Start_angle,
					     original_parameters->Opening_angle,
					     original_parameters->numSide,
					     original_parameters->Num_z_planes,
					     original_parameters->Z_values,
					     original_parameters->Rmin,
					     original_parameters->Rmax);
	}
	else {
		G4cerr << "G4Polyhedra: visualization of this type of G4Polyhedra is not supported at this time" << endl;
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




//
// G4Polyhedra::G4PolyhedraHistorical stuff
//
G4Polyhedra::G4PolyhedraHistorical::~G4PolyhedraHistorical()
{
	delete [] Z_values;
	delete [] Rmin;
	delete [] Rmax;
}

G4Polyhedra::G4PolyhedraHistorical::G4PolyhedraHistorical( const G4PolyhedraHistorical &source )
{
	Start_angle 	= source.Start_angle;
	Opening_angle	= source.Opening_angle;
	numSide		= source.numSide;
	Num_z_planes	= source.Num_z_planes;
	
	Z_values	= new G4double[Num_z_planes];
	Rmin		= new G4double[Num_z_planes];
	Rmax		= new G4double[Num_z_planes];
	
	G4int i;
	for( i = 0; i < Num_z_planes; i++) {
		Z_values[i] = source.Z_values[i];
		Rmin[i]	    = source.Rmin[i];
		Rmax[i]	    = source.Rmax[i];
	}
}

