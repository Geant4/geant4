// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polycone.cc,v 1.1 2000-04-07 11:01:50 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4Polycone.cc
//
// Implementation of a CSG polycone
//
// --------------------------------------------------------------------

#include "G4Polycone.hh"
#include "G4PolyconeSide.hh"
#include "G4PolyPhiFace.hh"

#include "G4Polyhedron.hh"
#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"


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
	// Some historical ugliness
	//
	original_parameters = new G4PolyconeHistorical();
	
	original_parameters->Start_angle = phiStart;
	original_parameters->Opening_angle = phiTotal;
	original_parameters->Num_z_planes = numZPlanes;
	original_parameters->Z_values = new G4double[numZPlanes];
	original_parameters->Rmin = new G4double[numZPlanes];
	original_parameters->Rmax = new G4double[numZPlanes];
        G4int i;
	for (i=0; i<numZPlanes; i++) {
		original_parameters->Z_values[i] = zPlane[i];
		original_parameters->Rmin[i] = rInner[i];
		original_parameters->Rmax[i] = rOuter[i];
	}
		
	//
	// Build RZ polygon using special PCON/PGON GEANT3 constructor
	//
	G4ReduciblePolygon *rz = new G4ReduciblePolygon( rInner, rOuter, zPlane, numZPlanes );
	
	//
	// Do the real work
	//
	Create( phiStart, phiTotal, rz );
	
	delete rz;
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
	original_parameters = 0;

	G4ReduciblePolygon *rz = new G4ReduciblePolygon( r, z, numRZ );
	
	Create( phiStart, phiTotal, rz );
	
	delete rz;
}


//
// Create
//
// Generic create routine, called by each constructor after conversion of arguments
//
void G4Polycone::Create( const G4double phiStart,
            	     	 const G4double phiTotal,
		     	 G4ReduciblePolygon *rz	  )
{
	//
	// Perform checks of rz values
	//
	if (rz->Amin() < 0.0) 
		G4Exception( "G4Polycone: Illegal input parameters: All R values must be >= 0" );
		
	G4double rzArea = rz->Area();
	if (rzArea < -kCarTolerance) rz->ReverseOrder();

	else if (rzArea < -kCarTolerance)
		G4Exception( "G4Polycone: Illegal input parameters: R/Z cross section is zero or near zero" );
		
	if ((!rz->RemoveDuplicateVertices( kCarTolerance )) || 
	    (!rz->RemoveRedundantVertices( kCarTolerance ))     ) 
	    	G4Exception( "G4Polycone: Illegal input parameters: Too few unique R/Z values" );

	if (rz->CrossesItself(1/kInfinity)) 
		G4Exception( "G4Polycone: Illegal input parameters: R/Z segments cross" );

	numCorner = rz->NumVertices();

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
	// Allocate corner array. 
	//
	corners = new G4PolyconeSideRZ[numCorner];

	//
	// Copy corners
	//
	G4ReduciblePolygonIterator iterRZ(rz);
	
	G4PolyconeSideRZ *next = corners;
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
		
		//
		// We must decide here if we can dare declare one of our faces
		// as having a "valid" normal (i.e. allBehind = true). This
		// is never possible if the face faces "inward" in r.
		//
		G4bool allBehind;
		if (corner->z > next->z) {
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
		
		*face++ = new G4PolyconeSide( prev, corner, next, nextNext,
					      startPhi, endPhi-startPhi, phiIsOpen, allBehind );
	} while( prev=corner, corner=next, corner > corners );
	
	if (phiIsOpen) {
		//
		// Construct phi open edges
		//
		*face++ = new G4PolyPhiFace( rz, startPhi, 0, endPhi  );
		*face++ = new G4PolyPhiFace( rz, endPhi,   0, startPhi );
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
G4Polycone::~G4Polycone()
{
	delete [] corners;
	
	if (original_parameters) delete original_parameters;
}




//
// Copy constructor
//
G4Polycone::G4Polycone( const G4Polycone &source ) : G4VCSGfaceted( source )
{
	CopyStuff( source );
}


//
// Assignment operator
//
const G4Polycone &G4Polycone::operator=( const G4Polycone &source )
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
void G4Polycone::CopyStuff( const G4Polycone &source )
{
	//
	// Simple stuff
	//
	startPhi	= source.startPhi;
	endPhi		= source.endPhi;
	phiIsOpen	= source.phiIsOpen;
	numCorner	= source.numCorner;

	//
	// The corner array
	//
	corners = new G4PolyconeSideRZ[numCorner];
	
	G4PolyconeSideRZ	*corn = corners,
				*sourceCorn = source.corners;
	do {
		*corn = *sourceCorn;
	} while( ++sourceCorn, ++corn < corners+numCorner );
	
	//
	// Original parameters
	//
	if (source.original_parameters) {
		original_parameters = new G4PolyconeHistorical( *source.original_parameters );
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
EInside G4Polycone::Inside( const G4ThreeVector &p ) const
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
G4double G4Polycone::DistanceToIn( const G4ThreeVector &p, const G4ThreeVector &v ) const
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
	// This has to be fixed in visualization. Fake it for the moment.
	// 
	if (original_parameters) {
	
		return new G4PolyhedronPcon( original_parameters->Start_angle,
					     original_parameters->Opening_angle,
					     original_parameters->Num_z_planes,
					     original_parameters->Z_values,
					     original_parameters->Rmin,
					     original_parameters->Rmax);
	}
	else {
		G4cerr << "G4Polycone: visualization of this type of G4Polycone is not supported at this time" << G4endl;
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


//
// G4Polycone:G4PolyconeHistorical stuff
//
G4Polycone::G4PolyconeHistorical::~G4PolyconeHistorical()
{
	delete [] Z_values;
	delete [] Rmin;
	delete [] Rmax;
}

G4Polycone::G4PolyconeHistorical::G4PolyconeHistorical( const G4Polycone::G4PolyconeHistorical &source )
{
	Start_angle 	= source.Start_angle;
	Opening_angle	= source.Opening_angle;
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

