//
// FredTest3
//
// Implementation of fred's third geometry tester
//
#include "FredTest3.hh"

#include "Randomize.hh"
#include "G4VSolid.hh"
#include "G4UImanager.hh"

#include <time.h>
#include "g4std/iomanip"
#include "g4std/strstream"

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
		 
		 
	if (grids.x() > 0) dx = grids.x()*rint( dx/grids.x() );
	if (grids.y() > 0) dy = grids.y()*rint( dx/grids.y() );
	if (grids.z() > 0) dz = grids.z()*rint( dx/grids.z() );

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
	if (cutoff <= 0) G4Exception( "Illegal cutoff" );

	G4double answer;
	do {
		answer = -3.0;
		for( G4int j = 0; j < 6; j++ ) answer += G4UniformRand();
		answer *= sqrt(2.0);
	} while( fabs(answer) > cutoff );
	
	return(answer);
}


//
// RunTest
//
// Do your stuff!
//
void FredTest3::RunTest( const G4VSolid *testVolume, G4std::ostream &logger )
{
	//
	// Clear error list
	//
	ClearErrors();
	
	//
	// Output test parameters
	//
	time_t now = time(0);
        char timebuf[25];
        timebuf[24]=0;
        strncpy( ctime(&now), timebuf, 24 );
	G4String dateTime( timebuf );
	
	logger << "% Fred test3 logged output " << dateTime << G4endl;
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
	logger << G4std::setprecision(14);
	
	//
	// Set clock
	//
	clock_t start = clock();
	
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

	now = time(0);
        char timebuf2[25];
        timebuf2[24]=0;
        strncpy( ctime(&now), timebuf2, 24 );
	G4String dateTime2( timebuf2 );
	logger << dateTime2 << G4endl;

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
G4int FredTest3::DebugError( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
{
	G4ThreeVector p, v;
	
	//
	// Recover information from log file
	//
	G4int error = GetLoggedPV( logger, errorIndex, p, v );
	if (error) return error;
			
	//
	// Setup fred to simulate this event
	//
	// We do this using the command line, a cheesy short-cut but
	// rather useful in hacks like this.
	//
	// If you are writing your own serious GEANT4 application,
	// please do something better.
	//
	char commandBuffer[255];
	G4UImanager *UI = G4UImanager::GetUIpointer();

	UI->ApplyCommand( "/fred/gun G4" );
	
	UI->ApplyCommand( "/gun/particle geantino" );
	
	G4std::ostrstream formatter1( commandBuffer, 255 );
	formatter1 << G4std::setprecision(14);
	formatter1 << "/gun/position "  << p.x() << " " << p.y() << " " << p.z() << " mm" << G4endl;
	UI->ApplyCommand( commandBuffer );
	
	G4std::ostrstream formatter2( commandBuffer, 255 );
	formatter2 << G4std::setprecision(14);
	formatter2 << "/gun/direction " << v.x() << " " << v.y() << " " << v.z() << " mm" << G4endl;
	UI->ApplyCommand( commandBuffer );
	return 0;
}


//
// DebugInside
//
// Recover previously logged error and invoke G4VSolid::Inside
//
G4int FredTest3::DebugInside( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
	EInside answer = testVolume->Inside( p );
	return 0;
}


//
// DebugToInP
//
// Recover previously logged error and invoke G4VSolid::DistanceToIn(p)
//
G4int FredTest3::DebugToInP( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
	G4double answer = testVolume->DistanceToIn( p );
	return 0;
}


//
// DebugToInPV
//
// Recover previously logged error and invoke G4VSolid::DistanceToIn(p,v)
//
G4int FredTest3::DebugToInPV( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
	
	p += answer*v;
	
	EInside inside = testVolume->Inside(p);
	return 0;
}


//
// DebugToOutP
//
// Recover previously logged error and invoke G4VSolid::DistanceToOut(p)
//
G4int FredTest3::DebugToOutP( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
	G4double answer = testVolume->DistanceToOut( p );
	return 0;
}


//
// DebugToOutPV
//
// Recover previously logged error and invoke G4VSolid::DistanceToOut(p,v)
//
G4int FredTest3::DebugToOutPV( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
	G4double answer = testVolume->DistanceToOut( p, v );
	
	p += answer*v;
	
	EInside inside = testVolume->Inside(p);
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
				  const FredTest3PointList *inside, const G4ThreeVector point, G4std::ostream &logger )
{
	G4int i, n = inside->NumPoints();
	
	G4double safeDistance = testVolume->DistanceToIn( point );
	if (safeDistance <= 0.0) {
		ReportError( nError, point, 0, "T0: DistanceToIn(p) <= 0", logger );
		return;
	}
	
	for( i=0; i < n; i++ ) {
		G4ThreeVector vr = (*inside)[i] - point;
		G4ThreeVector v = vr.unit();
		
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
		
		G4ThreeVector p = point + dist*v;
		
		EInside insideOrNot = testVolume->Inside( p );
		if (insideOrNot == kOutside) {
			ReportError( nError, point, v, "T0: DistanceToIn(p,v) undershoots", logger );
			continue;
		}
		if (insideOrNot == kInside) {
			ReportError( nError, point, v, "TO: DistanceToIn(p,v) overshoots", logger );
			continue;
		}
		
		dist = testVolume->DistanceToIn( p );
		if (dist != 0) {
			ReportError( nError, p, v, "T02: DistanceToIn(p) should be zero", logger );
			continue;
		}
		
		dist = testVolume->DistanceToOut( p );
		if (dist != 0) {
			ReportError( nError, p, v, "T02: DistanceToOut(p) should be zero", logger );
			continue;
		}
		
		dist = testVolume->DistanceToIn( p, v );
		if (dist != 0) {
			ReportError( nError, p, v, "T02: DistanceToIn(p,v) should be zero", logger );
			continue;
		}	
		
		G4bool validNorm;
		G4ThreeVector norm;
		
		dist = testVolume->DistanceToOut( p, v, true, &validNorm, &norm );
		if (dist == 0) continue;
		
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
				 const FredTest3PointList *outside, const G4ThreeVector point, G4std::ostream &logger )
{
	G4int i, n = outside->NumPoints();
	
	G4double safeDistance = testVolume->DistanceToOut( point );
	if (safeDistance <= 0.0) {
		ReportError( nError, point, 0, "TI: DistanceToOut(p) <= 0", logger );
		return;
	}
	
	for( i=0; i < n; i++ ) {
		G4ThreeVector vr = (*outside)[i] - point;
		G4ThreeVector v = vr.unit();

		G4bool validNorm;
		G4ThreeVector norm;
		
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
			ReportError( nError, point, v, "TI: DistanceToOut(p,v) < DistanceToIn(p)", logger );
			continue;
		}

		if (validNorm) {
			if (norm.dot(v) < 0) {
				ReportError( nError, point, v, "TI: Outgoing normal incorrect", logger );
				continue;
			}
		}
		
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
			     const G4ThreeVector v, const G4String comment, G4std::ostream &logger )
{
	//
	// Have we encountered this error message before?
	//
	FredTest3ErrorList *last, *errors = errorList;
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
G4int FredTest3::GetLoggedPV( G4std::istream &logger, const G4int errorIndex,
			      G4ThreeVector &p, G4ThreeVector &v        ) const
{
	logger >> G4std::setprecision(14);		// I wonder if this is necessary?

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
	delete pointList;
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
		G4int	irand = G4UniformRand()*( (G4double)maxPoints );
		pointList[irand] = newPoint;
	}
}
