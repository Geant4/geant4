//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// FredTest3
//
// Implementation of fred's third geometry tester
//
#include "FredTest3.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include "G4VSolid.hh"
#include "G4UImanager.hh"
#include "G4GeometryTolerance.hh"

#include <time.h>
#include <iomanip>
#include <sstream>

//
// Constructor
//
FredTest3::FredTest3()
{
	//
	// Defaults
	//
	SetDefaults();
	
	//
	// Zero error list
	//
	errorList = 0;
}


//
// Destructor
//
FredTest3::~FredTest3()
{
	ClearErrors();
}


//
// SetDefaults
//
// Set defaults values
//
void FredTest3::SetDefaults()
{
	target = G4ThreeVector( 0, 0, 0 );
	widths = G4ThreeVector( 1*m, 1*m, 1*m );
	grids  = G4ThreeVector( 0, 0, 0 );
	
	maxPoints = 10000;
	maxErrors = 100;
}


//
// GetRandomPoint
//
// Return a random point in three dimensions using the current
// specs 
//
G4ThreeVector FredTest3::GetRandomPoint() const {
	G4double dx = widths.x()*GaussianRandom(10*m/widths.x()),
		 dy = widths.y()*GaussianRandom(10*m/widths.y()),
		 dz = widths.z()*GaussianRandom(10*m/widths.z());
		 
		 
	if (grids.x() > 0) dx = grids.x()*G4int( dx/grids.x()+1 );
	if (grids.y() > 0) dy = grids.y()*G4int( dx/grids.y()+1 );
	if (grids.z() > 0) dz = grids.z()*G4int( dx/grids.z()+1 );

	G4ThreeVector randvec( dx, dy, dz );
	
	return target + randvec;
}


//
// GaussianRandom
//
// Return Gaussian random number of unit width
//
// A classic, slow, but remarkably effective algorithm. Certainly good
// enough for our purposes.
//
G4double FredTest3::GaussianRandom(const G4double cutoff) const {
	if (cutoff <= 0) 
	  G4Exception("FredTest3","Fred001",FatalException,"Illegal cutoff" );

	G4double answer;
	do {
		answer = -3.0;
		for( G4int j = 0; j < 6; j++ ) answer += G4UniformRand();
		answer *= std::sqrt(2.0);
	} while( std::fabs(answer) > cutoff );
	
	return(answer);
}


//
// RunTest
//
// Do your stuff!
//
void FredTest3::RunTest( const G4VSolid *testVolume, std::ostream &logger )
{
	//
	// Clear error list
	//
	ClearErrors();
	
	//
	// Output test parameters
	//
	time_t now;
	time(&now);
	G4String dateTime(ctime(&now));
	
	logger << "% Fred test3 logged output " << dateTime;
	logger << "% volume =    " << testVolume->GetName() << G4endl;
	logger << "% target =    " << target << G4endl;
	logger << "% widths =    " << widths << G4endl;
	logger << "% grids  =    " << grids  << G4endl;
	logger << "% maxPoints = " << maxPoints << G4endl;
	logger << "% maxErrors = " << maxErrors << G4endl;

	//
	// Setup lists
	//
	FredTest3PointList inside(100);
	FredTest3PointList outside(100);
	FredTest3PointList surface(100);
	
	//
	// Set iostream precision to 14 digits
	//
	logger << std::setprecision(14);
	
	//
	// Set clock
	//
	// clock_t start = clock();
	
	//
	// Loop over points
	//
	G4int nIn = 0, nOut = 0, nSurf = 0;
	
	G4int nPoint = 0;
	G4int nError = 0;
	
	for(;;) {
  
		//
		// Generate a random point
		//
		G4ThreeVector	point = GetRandomPoint();
		
		//
		// Catagorize, test, and store
		//
		EInside catagory = testVolume->Inside( point );
		switch( catagory ) {
			case kOutside:
			nOut++;
			TestOutsidePoint( testVolume, &nError, &inside, point, logger );
			outside.AddPoint( point );
			break;
			
			case kInside:
			nIn++;
			TestInsidePoint( testVolume, &nError, &outside, point, logger );
			inside.AddPoint( point );
			break;
			
			case kSurface:
			nSurf++;
//			TestSurfacePoint( testVolume, &nError, point, logger );
			surface.AddPoint( point );
			break;
		}
		
		//
		// Return early if we have accumulated enough errors
		//
		if (nError >= maxErrors) {
			logger << "% End of test (maximum number errors) ";
			break;
		}
		if (++nPoint >= maxPoints) {
			logger << "% End of test (maximum number points) ";
			break;
		}
	
	}

	time(&now);
	G4String dateTime2(ctime(&now));
	logger << dateTime2;

	logger << "% Statistics: points=" << nPoint << " errors reported=" << nError << G4endl;

	logger << "%             inside=" << nIn << " outside=" << nOut
               << " surface=" << nSurf << G4endl;
        logger << "%             cpu time=" << clock()/CLOCKS_PER_SEC << G4endl;
	       
	logger << "%(End of file)" << G4endl;
}


//
// DebugError
//
// Recover previously logged error and setup particle gun appropriately
//
void FredTest3::RunDebug( const G4VSolid* solid, std::istream &logger )
{

	G4UImanager *UI = G4UImanager::GetUIpointer();
        G4int errorIndex = 0;
        while ( ! DebugError(solid, logger, ++errorIndex) ) {

 	        UI->ApplyCommand( "/run/beamOn 1" );
        }        
}

//
// DebugError
//
// Recover previously logged error and setup particle gun appropriately
//
G4int FredTest3::DebugError( const G4VSolid *, std::istream &logger,
                             const G4int errorIndex ) const
{
	G4ThreeVector p, v;
	
	//
	// Recover information from log file
	//
	G4int error = GetLoggedPV( logger, errorIndex, p, v );
	if (error) return error;
	
        G4cout << "DebugError  " << errorIndex << ":  p=" << p << ",  v=" << v << G4endl; 
			
	//
	// Setup fred to simulate this event
	//
	// We do this using the command line, a cheesy short-cut but
	// rather useful in hacks like this.
	//
	// If you are writing your own serious GEANT4 application,
	// please do something better.
	//
	G4UImanager *UI = G4UImanager::GetUIpointer();

	UI->ApplyCommand( "/tracking/verbose 1" );

	UI->ApplyCommand( "/fred/gun G4" );
	
	UI->ApplyCommand( "/gun/particle geantino" );
	
	std::ostringstream formatter1;
	formatter1 << std::setprecision(14);
	formatter1 << "/gun/position "  << p.x() << " " << p.y() << " " << p.z() << " mm" << G4endl;
        G4String commandBuffer=formatter1.str();
	UI->ApplyCommand( commandBuffer );
	
	std::ostringstream formatter2;
	formatter2 << std::setprecision(14);
	formatter2 << "/gun/direction " << v.x() << " " << v.y() << " " << v.z() << " mm" << G4endl;
        commandBuffer=formatter2.str();
	UI->ApplyCommand( commandBuffer );
	return 0;
}


//
// DebugInside
//
// Recover previously logged error and invoke G4VSolid::Inside
//
G4int FredTest3::DebugInside( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const
{	
	G4ThreeVector p, v;
	
	//
	// Recover information from log file
	//
	G4int error = GetLoggedPV( logger, errorIndex, p, v );
	if (error) return error;
	
	//
	// Call
	//
	G4cout << "testVolume->Inside(p): " << testVolume->Inside( p ) << G4endl; 
	return 0;
}


//
// DebugToInP
//
// Recover previously logged error and invoke G4VSolid::DistanceToIn(p)
//
G4int FredTest3::DebugToInP( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const
{	
	G4ThreeVector p, v;
	
	//
	// Recover information from log file
	//
	G4int error = GetLoggedPV( logger, errorIndex, p, v );
	if (error) return error;
	
	//
	// Call
	//
	G4cout << "testVolume->DistanceToIn(p): " <<  testVolume->DistanceToIn( p )  << G4endl; 
	return 0;
}


//
// DebugToInPV
//
// Recover previously logged error and invoke G4VSolid::DistanceToIn(p,v)
//
G4int FredTest3::DebugToInPV( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const
{	
	G4ThreeVector p, v;
	
	//
	// Recover information from log file
	//
	G4int error = GetLoggedPV( logger, errorIndex, p, v );
	if (error) return error;
	
	//
	// Call
	//
	G4double answer = testVolume->DistanceToIn( p, v );
	G4cout << "testVolume->DistanceToIn(p,v): " << answer << G4endl; 
	
	p += answer*v;
	
	G4cout << "testVolume->Inside(p+=answer*v):" <<  p << " " << testVolume->Inside(p) << G4endl;
	return 0;
}


//
// DebugToOutP
//
// Recover previously logged error and invoke G4VSolid::DistanceToOut(p)
//
G4int FredTest3::DebugToOutP( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const
{	
	G4ThreeVector p, v;
	
	//
	// Recover information from log file
	//
	G4int error = GetLoggedPV( logger, errorIndex, p, v );
	if (error) return error;
	
	//
	// Call
	//
	G4cout << "testVolume->DistanceToOut(p): " << testVolume->DistanceToOut( p ) << G4endl; 

	return 0;
}


//
// DebugToOutPV
//
// Recover previously logged error and invoke G4VSolid::DistanceToOut(p,v)
//
G4int FredTest3::DebugToOutPV( const G4VSolid *testVolume, std::istream &logger, const G4int errorIndex ) const
{	
	G4ThreeVector p, v;
	
	//
	// Recover information from log file
	//
	G4int error = GetLoggedPV( logger, errorIndex, p, v );
	if (error) return error;

	//
	// Call
	//
	G4bool validNorm;
	G4ThreeVector norm;
	G4double answer = testVolume->DistanceToOut( p, v, true, &validNorm, &norm);
	G4cout << "testVolume->DistanceToOut( p, v ): " << answer << " validNorm: " << validNorm << G4endl; 
	
	p += answer*v;
	G4cout << "testVolume->Inside(p += answer*v): " << testVolume->Inside(p) << G4endl; 
	
	return 0;
}


//
// --------------------------------------
// Tests
//

//
// TestOutsidePoint
//
void FredTest3::TestOutsidePoint( const G4VSolid *testVolume, G4int *nError,
				  const FredTest3PointList *inside, const G4ThreeVector point, std::ostream &logger )
{
	// Check if point is really outside
	G4double safeDistance = testVolume->DistanceToIn( point );
	if (safeDistance <= 0.0) {
		ReportError( nError, point, G4ThreeVector(), "T0: DistanceToIn(p) <= 0", logger );
		return;
	}
	
	// Loop over inside points
	for( G4int i=0; i < inside->NumPoints(); i++ ) {
  
  		// Direction vector from tested point to the inside point
		G4ThreeVector vr = (*inside)[i] - point;
		G4ThreeVector v = vr.unit();
		
  		// Distance from the outside point to solid in direction
  		// to the tested inside point
		G4double dist = testVolume->DistanceToIn( point, v );

		if (dist <= 0) {
			ReportError( nError, point, v, "T0: DistanceToIn(p,v) <= 0", logger );
			continue;
		}
		if (dist == kInfinity) {
			ReportError( nError, point, v, "T0: DistanceToIn(p,v) == kInfinity", logger );
			continue;
		}
		if (dist < safeDistance+1E-10) {
			ReportError( nError, point, v, "T0: DistanceToIn(p,v) < DistanceToIn(p)", logger );
			continue;
		}
		
  		// Tested point moved to the solid surface
  		// (in direction to the tested inside point): 
		G4ThreeVector p = point + dist*v;
		
  		// Inside(p) should return kSurface.
  		// If kOutside is returned, report undershoots; if kInside, report overshoots
		EInside insideOrNot = testVolume->Inside( p );
		if (insideOrNot == kOutside) {
			ReportError( nError, point, v, "T0: DistanceToIn(p,v) undershoots", logger );
			continue;
		}
		if (insideOrNot == kInside) {
			ReportError( nError, point, v, "TO: DistanceToIn(p,v) overshoots", logger );
			continue;
		}
		
  		// DistanceToIn for the computed point on the surface
		dist = testVolume->DistanceToIn( p );
		//if (dist != 0) {
		if (dist > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {
			ReportError( nError, p, v, "T02: DistanceToIn(p) should be zero", logger );
			continue;
		}
		
  		// DistanceToOut for the computed point on the surface
		dist = testVolume->DistanceToOut( p );
		//if (dist != 0) {
		if (dist > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {
			ReportError( nError, p, v, "T02: DistanceToOut(p) should be zero", logger );
			continue;
		}
		
  		// DistanceToIn with direction for the computed point on the surface
		dist = testVolume->DistanceToIn( p, v );
		//if (dist != 0) {
		if (dist > G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) {
			ReportError( nError, p, v, "T02: DistanceToIn(p,v) should be zero", logger );
			continue;
		}	
		
  		// DistanceToOut with direction for the computed point on the surface
		G4bool validNorm;
		G4ThreeVector norm;
		
		dist = testVolume->DistanceToOut( p, v, true, &validNorm, &norm );
		// if (dist == 0) continue;
		// if (dist < G4GeometryTolerance::GetInstance()->GetSurfaceTolerance()) continue;
       		// Why quit here withour Error ??? 
    
		if (dist == kInfinity) {
			ReportError( nError, p, v, "T02: DistanceToOut(p,v) == kInfinity", logger );
			continue;
		}
		
		if (validNorm) {
			if (norm.dot(v) < 0) {
				ReportError( nError, p, v, "T02: Outgoing normal incorrect", logger );
				continue;
			}
		}
		
  		// The point on surface (entering point) moved to the second point
  		// on the surface (exiting point)
  		// (in direction to the tested inside point):
		G4ThreeVector p2 = p + v*dist;
		
		insideOrNot = testVolume->Inside(p2);
		if (insideOrNot == kInside) {
			ReportError( nError, p, v, "T02: DistanceToOut(p,v) undershoots", logger );
			continue;
		}
		if (insideOrNot == kOutside) {
			ReportError( nError, p, v, "TO2: DistanceToOut(p,v) overshoots", logger );
			continue;
		}
			
  		// DistanceToIn from the exiting point (in the exiting direction)
  		// If infinity, there is no intersectin with solid anymore, what means that the solid
  		// lied entirely between entering and existing points, and the validNorm should have been true.
                
		dist = testVolume->DistanceToIn(p2,v);
		if (validNorm) {
			if (dist != kInfinity) {
				ReportError( nError, p, v, "TO2: DistanceToOut incorrectly returns validNorm==true", logger );
				continue;
			}
		}
		else {
			if (dist == kInfinity) {
				ReportError( nError, p, v, "TO2: DistanceToOut incorrectly returns validNorm==false", logger );
				continue;
			}
		}
	}
}



//
// TestInsidePoint
//
void FredTest3::TestInsidePoint( const G4VSolid *testVolume, G4int *nError,
				 const FredTest3PointList *outside, const G4ThreeVector point, std::ostream &logger )
{
	// Check if point is really inside
	G4double safeDistance = testVolume->DistanceToOut( point );
	if (safeDistance <= 0.0) {
		ReportError( nError, point, G4ThreeVector(), "TI: DistanceToOut(p) <= 0", logger );
		return;
	}
	
	// Loop over outside points
	for( G4int i=0; i < outside->NumPoints(); i++ ) {

  		// Direction vector from tested point to the outside point
		G4ThreeVector vr = (*outside)[i] - point;
		G4ThreeVector v = vr.unit();

		G4bool validNorm;
		G4ThreeVector norm;
		
  		// Distance from the inside point to solid border in direction
  		// to the tested outside point
		G4double dist = testVolume->DistanceToOut( point, v, true, &validNorm, &norm );
		if (dist <= 0) {
			ReportError( nError, point, v, "TI: DistanceToOut(p,v) <= 0", logger );
			continue;
		}
		if (dist == kInfinity) {
			ReportError( nError, point, v, "TI: DistanceToOut(p,v) == kInfinity", logger );
			continue;
		}
		if (dist < safeDistance+1E-10) {
			ReportError( nError, point, v, "TI: DistanceToOut(p,v) < DistanceToOut(p)", logger );
			continue;
		}

		if (validNorm) {
			if (norm.dot(v) < 0) {
				ReportError( nError, point, v, "TI: Outgoing normal incorrect", logger );
				continue;
			}
		}
		
  		// Tested point moved to the solid surface
  		// (in direction to the tested outside point): 
		G4ThreeVector p = point + v*dist;
		
		EInside insideOrNot = testVolume->Inside(p);
		if (insideOrNot == kInside) {
			ReportError( nError, point, v, "TI: DistanceToOut(p,v) undershoots", logger );
			continue;
		}
		if (insideOrNot == kOutside) {
			ReportError( nError, point, v, "TI: DistanceToOut(p,v) overshoots", logger );
			continue;
		}
	}
}


//
// ReportError
//
// Report the specified error message, but only if it has not been reported a zillion
// times already.
//
void FredTest3::ReportError( G4int *nError, const G4ThreeVector p, 
			     const G4ThreeVector v, const G4String comment, std::ostream &logger )
{
	//
	// Have we encountered this error message before?
	//
	FredTest3ErrorList *last=0, *errors = errorList;
	while( errors ) {
		//
		// Note: below is an expensive comparison. This could be replaced with something
		// faster if we really wanted, like a hash. But is it worth the trouble? Nah.
		//
		if (errors->message == comment) {
			//
			// Yup. Increment count, and return now if the count is too high
			//
			if ( ++errors->nUsed > 5 ) return;
			break;
		}
		last = errors;
		errors = errors->next;
	}
	
	if (errors == 0) {
		//
		// New error: add it the end of our list
		//
		errors = new FredTest3ErrorList;
		errors->message = comment;
		errors->nUsed = 1;
		errors->next = 0;
		if (errorList) 
			last->next = errors;
		else
			errorList = errors;
	}

	//
	// Output the message
	//	
	logger << "% " << comment;
	if (errors->nUsed == 5) logger << " (any further such errors suppressed)";
	logger << G4endl;
	
	logger << ++(*nError) << " " << p.x() << " " << p.y() << " " << p.z() 
			      << " " << v.x() << " " << v.y() << " " << v.z() << G4endl;
}


//
// ClearErrors
//
// Reset list of errors (and clear memory)
//
void FredTest3::ClearErrors()
{
	FredTest3ErrorList *here, *next;
	
	here = errorList;
	while( here ) {
		next = here->next;
		delete here;
		here = next;
	}
}


//
// GetLoggedPV
//
// Get the p and v vectors stored in a test3 log file
//
G4int FredTest3::GetLoggedPV( std::istream &logger, const G4int errorIndex,
			      G4ThreeVector &p, G4ThreeVector &v        ) const
{
	logger >> std::setprecision(14);		// I wonder if this is necessary?
        logger.seekg(0);

	//
	// Search for the requested error index, skipping comments along the way
	//
	for(;;) {
		while( logger.peek() == '%' ) logger.ignore( 99999, '\n' );
		
		if (logger.peek() == EOF) return 1;
		
		G4int thisIndex;
		
		logger >> thisIndex;
                
		if (thisIndex == errorIndex) break;
		
		logger.ignore( 99999, '\n' );
	}
	
	//
	// Extract the vectors
	//
	// We should probably have some error checking here...
	//
	G4double x, y, z;
	
	logger >> x >> y >> z;
	p = G4ThreeVector(x*mm,y*mm,z*mm);
	
	logger >> x >> y >> z;
	v = G4ThreeVector(x*mm,y*mm,z*mm);

	return 0;
}


//
// --------------------------------------
// FredTest3PointList stuff
//

FredTest3PointList::FredTest3PointList( G4int size )
{
	pointList = new G4ThreeVector[size];
	maxPoints = size;
	numPoints = 0;
}

FredTest3PointList::~FredTest3PointList()
{
	delete [] pointList;
}

void FredTest3PointList::AddPoint( G4ThreeVector newPoint )
{
	if (numPoints < maxPoints) {
		//
		// Since we still have space, append the point on the
		// end of the list
		//
		pointList[numPoints] = newPoint;
		numPoints++;
	}
	else {
		//
		// Our list is filled up, so replace a random
		// entry
		//
		G4int irand = G4int(G4UniformRand()*( (G4double)maxPoints ));
		pointList[irand] = newPoint;
	}
}
