//
// SBTrun
//
// Implementation of the batch solid test
//
#include "SBTrun.hh"

#include "Randomize.hh"
#include "G4VSolid.hh"

#include "SBTVisManager.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"
#include "g4std/iomanip"
#include "g4std/strstream"

#include "G4SolidStore.hh"

#include <time.h>

#define DEBUG 1

//
// Constructor
//
SBTrun::SBTrun()
{
  //
  // Defaults
  //
  SetDefaults();

  //
  // Zero error list
  //
  errorList = 0;
  CurrentSolid = "";
}


//
// Destructor
//
SBTrun::~SBTrun()
{
  ClearErrors();
}


//
// SetDefaults
//
// Set defaults values
//
void SBTrun::SetDefaults()
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
G4ThreeVector SBTrun::GetRandomPoint() const {
  G4double dx = widths.x()*GaussianRandom(10*m/widths.x()),
    dy = widths.y()*GaussianRandom(10*m/widths.y()),
    dz = widths.z()*GaussianRandom(10*m/widths.z());


  if (grids.x() > 0) dx = grids.x()*rint( dx/grids.x() );
  if (grids.y() > 0) dy = grids.y()*rint( dy/grids.y() );
  if (grids.z() > 0) dz = grids.z()*rint( dz/grids.z() );

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
G4double SBTrun::GaussianRandom(const G4double cutoff) const {
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
void SBTrun::RunTest( const G4VSolid *testVolume, G4std::ostream &logger )
{
  //
  // Clear error list
  //
    
  ClearErrors();

	  //
	  // Output test parameters
	  //
  time_t now = time(0);
  //	G4String dateTime( ctime(&now), 24 );		// AFAIK, this is standard c++
  //      it looks like the above is missing from STLInterface
  char timebuf[25];
  timebuf[24]=0;
  strncpy( ctime(&now), timebuf, 24 );
  G4String dateTime( timebuf );

  {
    /* 
       MEDERNACH Emmanuel
       aug 2000
	       
       We try to retrieve the solid being tested 
       to write in the log file
    */

    G4SolidStore* SolidStore = G4SolidStore::GetInstance();
    G4cout << (SolidStore->at(0))->GetName() << G4endl ;
  }


  logger << "% SBT logged output " << dateTime << G4endl;
  logger << "% " << CurrentSolid << G4endl;
  logger << "% target =    " << target << G4endl;
  logger << "% widths =    " << widths << G4endl;
  logger << "% grids  =    " << grids  << G4endl;
  logger << "% maxPoints = " << maxPoints << G4endl;
  logger << "% maxErrors = " << maxErrors << G4endl;

	  //
	  // Setup lists
	  //
  SBTrunPointList inside(100);
  SBTrunPointList outside(100);
  SBTrunPointList surface(100);

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

  G4std::setprecision(16);

  for(;;) {
    //
    // Generate a random point
    //
    G4ThreeVector	point = GetRandomPoint();

#if DEBUG
    G4cerr << G4std::setprecision(32);
    G4cerr << "Current point is : " << point* (1/cm) << G4endl ;
#endif
		  //
		  // Catagorize, test, and store
		  //
    EInside catagory = testVolume->Inside( point );


    switch( catagory ) {
    case kOutside:
      //G4cerr << "kOutside " << G4endl;
      
      nOut++;
      TestOutsidePoint( testVolume, &nError, &inside, &outside, point, logger );
      outside.AddPoint( point );
      break;

    case kInside:
      //G4cerr << "kInside " << G4endl;
      
      nIn++;
      TestInsidePoint( testVolume, &nError, &outside, point, logger );
      inside.AddPoint( point );
      break;

    case kSurface:
      //G4cerr << "kSurface " << G4endl;
      
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
  //	G4String dateTime2( ctime(&now), 24 );		
  strncpy( ctime(&now), timebuf, 24 );
  G4String dateTime2( timebuf );
  logger << dateTime2 << G4endl;

  logger << "% Statistics: points=" << nPoint << " errors=" << CountErrors()
	 << " errors reported=" << nError << G4endl;

  logger << "%             inside=" << nIn << " outside=" << nOut
	 << " surface=" << nSurf << G4endl;
  logger << "%             cpu time=" << clock()/CLOCKS_PER_SEC << G4endl;

  logger << "%(End of file)" << G4endl;
}


//
// DrawError
//
// Recover previously logged error and display it
//
G4int SBTrun::DrawError( const G4VSolid *testVolume, G4std::istream &logger, 
			 const G4int errorIndex, SBTVisManager *visManager ) const
{
  G4ThreeVector p, v;

  //
  // Recover information from log file
  //
  G4int error = GetLoggedPV( logger, errorIndex, p, v );
  if (error) return error;

  //
  // Draw away
  //
  visManager->ClearView();

  //
  // This draws the trajectory
  //    
  G4VisAttributes blueStuff( G4Color(0,0,1) );

  G4Polyline polyline;
  polyline.append( p );
  polyline.append( p + 4*v*m );
  polyline.SetVisAttributes( blueStuff );
  visManager->Draw( polyline );

  //
  // This draws the initial point p
  //
  G4Circle circle(p);
  circle.SetWorldSize( 5*cm );
  circle.SetVisAttributes( blueStuff );
  visManager->Draw( circle );

  //
  // This draws the target solid
  //
  G4VisAttributes redStuff( G4Color(1,0,0) );
  visManager->Draw( *testVolume, redStuff );

  visManager->Show();;

  return 0;
}


//
// DebugInside
//
// Recover previously logged error and invoke G4VSolid::Inside
//
G4int SBTrun::DebugInside( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
G4int SBTrun::DebugToInP( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
G4int SBTrun::DebugToInPV( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
G4int SBTrun::DebugToOutP( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
G4int SBTrun::DebugToOutPV( const G4VSolid *testVolume, G4std::istream &logger, const G4int errorIndex ) const
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
void SBTrun::TestOutsidePoint( const G4VSolid *testVolume, G4int *nError,
			       const SBTrunPointList *inside, const SBTrunPointList *outside,
			       const G4ThreeVector point, G4std::ostream &logger )
{
  G4int i, n = inside->NumPoints();
  
  G4double safeDistance = testVolume->DistanceToIn( point );

  if (safeDistance <= 0.0) {
    ReportError( nError, point, 0, safeDistance,"T0: DistanceToIn(p) <= 0", logger );
    return;
  }

  for( i=0; i < n; i++ ) {
    G4ThreeVector vr = (*inside)[i] - point;
    G4ThreeVector v = vr.unit();

    G4double dist = testVolume->DistanceToIn( point, v );
    if (dist <= 0) {
      ReportError( nError, point, v, safeDistance, "T0: DistanceToIn(p,v) <= 0", logger );
      continue;
    }
    if (dist == kInfinity) {
      ReportError( nError, point, v, safeDistance, "T0: DistanceToIn(p,v) == kInfinity", logger );
      continue;
    }
    if (dist < safeDistance-1E-10) {
      ReportError( nError, point, v, safeDistance, "T0: DistanceToIn(p,v) < DistanceToIn(p)", logger );
      continue;
    }

    G4ThreeVector p = point + dist*v;

    EInside insideOrNot = testVolume->Inside( p );
    if (insideOrNot == kOutside) {
      ReportError( nError, point, v, safeDistance, "T0: DistanceToIn(p,v) undershoots", logger );
      continue;
    }
    if (insideOrNot == kInside) {
      ReportError( nError, point, v, safeDistance, "TO: DistanceToIn(p,v) overshoots", logger );
      continue;
    }

    dist = testVolume->DistanceToIn( p );

    if (dist != 0) {
      ReportError( nError, p, v, safeDistance, "T02: DistanceToIn(p) should be zero", logger );
      // logger << "Dist != 0 : " << dist << endl;
      continue;
    }

    dist = testVolume->DistanceToOut( p );
    if (dist != 0) {
      ReportError( nError, p, v, safeDistance, "T02: DistanceToOut(p) should be zero", logger );
      continue;
    }

    dist = testVolume->DistanceToIn( p, v );
    safeDistance = testVolume->DistanceToIn( p );
    //
    // Beware! We might expect dist to be precisely zero, but this may not
    // be true at corners due to roundoff of the calculation of p = point + dist*v.
    // It should, however, *not* be infinity.
    //
    if (dist == kInfinity) {
      ReportError( nError, p, v, safeDistance, "T02: DistanceToIn(p,v) == kInfinity", logger );
      continue;
    }	

    G4bool validNorm;
    G4ThreeVector norm;

    dist = testVolume->DistanceToOut( p, v, true, &validNorm, &norm );
    if (dist == 0) continue;

    if (dist == kInfinity) {
      ReportError( nError, p, v, safeDistance, "T02: DistanceToOut(p,v) == kInfinity", logger );
      continue;
    }
    else if (dist < 0) {
      ReportError( nError, p, v, safeDistance, "T02: DistanceToOut(p,v) < 0", logger );
      continue;
    }

    if (validNorm) {
      if (norm.dot(v) < 0) {
	ReportError( nError, p, v, safeDistance, "T02: Outgoing normal incorrect", logger );
	continue;
      }
    }

    G4ThreeVector norm1 = testVolume->SurfaceNormal( p );
    if (norm1.dot(v) > 0) {
      ReportError( nError, p, v, safeDistance, "T02: Ingoing surfaceNormal is incorrect", logger );
    }


    G4ThreeVector p2 = p + v*dist;

    insideOrNot = testVolume->Inside(p2);
    if (insideOrNot == kInside) {
      ReportError( nError, p, v, safeDistance, "T02: DistanceToOut(p,v) undershoots", logger );
      continue;
    }
    if (insideOrNot == kOutside) {
      ReportError( nError, p, v, safeDistance, "TO2: DistanceToOut(p,v) overshoots", logger );
      continue;
    }

    G4ThreeVector norm2 = testVolume->SurfaceNormal( p2 );
    if (norm2.dot(v) < 0) {
      if (testVolume->DistanceToIn(p2,v) != 0)
	ReportError( nError, p2, v, safeDistance, "T02: Outgoing surfaceNormal is incorrect", logger );
    }
    if (validNorm) {
      if (norm.dot(norm2) < 0.0) {
	ReportError( nError, p2, v, safeDistance, "T02: SurfaceNormal and DistanceToOut disagree on normal", logger );
      }
    }

    if (validNorm) {
      dist = testVolume->DistanceToIn(p2,v);
      if (dist == 0) {
	//
	// We may have grazed a corner, which is a problem of design.
	// Check distance out
	//
	if (testVolume->DistanceToOut(p2,v) != 0) {
	  ReportError( nError, p, v, safeDistance, "TO2: DistanceToOut incorrectly returns validNorm==true (line of sight)(c)", logger );
	  continue;
	}
      }
      else if (dist != kInfinity) {
	ReportError( nError, p, v, safeDistance, "TO2: DistanceToOut incorrectly returns validNorm==true (line of sight)", logger );
	continue;
      }

      G4int j;
      for( j=0; j < n; j++ ) {
	G4ThreeVector p2top = (*inside)[j] - p2;

	if (p2top.dot(norm) > 0) {
	  ReportError( nError, p, v,safeDistance, 
		       "T02: DistanceToOut incorrectly returns validNorm==true (horizon)", logger );
	  continue;
	}
      }
    } // if valid normal
  } // Loop over inside points

  n = outside->NumPoints();

  for( i=0; i < n; i++ ) {
    G4ThreeVector vr = (*outside)[i] - point;
    if (vr.mag2() < DBL_MIN) continue;

    G4ThreeVector v = vr.unit();

#if DEBUG    
    G4cerr << G4std::setprecision(32);
    G4cerr << "point = " << point << G4endl ;
    G4cerr << "v = " << v << G4endl ;
#endif
    
    G4double dist = testVolume->DistanceToIn( point, v );

#if DEBUG
    G4cerr << "finished ..\n";
#endif    

    if (dist <= 0) {
      ReportError( nError, point, v, safeDistance, "T03: DistanceToIn(p,v) <= 0", logger );
      continue;
    }
    if (dist == kInfinity) {
      //G4cout << "dist == kInfinity" << G4endl ;
      continue;
    }
    if (dist < safeDistance-1E-10) {
      ReportError( nError, point, v, safeDistance, "T03: DistanceToIn(p,v) < DistanceToIn(p)", logger );
      continue;
    }
    G4ThreeVector p = point + dist*v;

    EInside insideOrNot = testVolume->Inside( p );
    if (insideOrNot == kOutside) {
      ReportError( nError, point, v, safeDistance, "T03: DistanceToIn(p,v) undershoots", logger );
      continue;
    }
    if (insideOrNot == kInside) {
      ReportError( nError, point, v, safeDistance, "TO3: DistanceToIn(p,v) overshoots", logger );
      continue;
    }
  } // Loop over outside points
}

//
// TestInsidePoint
//
void SBTrun::TestInsidePoint( const G4VSolid *testVolume, G4int *nError,
			      const SBTrunPointList *outside, const G4ThreeVector point, G4std::ostream &logger )
{
  G4int i, n = outside->NumPoints();

  G4double safeDistance = testVolume->DistanceToOut( point );
  if (safeDistance <= 0.0) {
    ReportError( nError, point, 0, safeDistance, "TI: DistanceToOut(p) <= 0", logger );
    return;
  }

  G4VisExtent visExtent(testVolume->GetExtent());
  if (point.x() < visExtent.GetXmin() ||
      point.x() > visExtent.GetXmax() || 
      point.y() < visExtent.GetYmin() || 
      point.y() > visExtent.GetYmax() || 
      point.z() < visExtent.GetZmin() || 
      point.z() > visExtent.GetZmax() ) {
    ReportError( nError, point, 0, safeDistance, "TI: Point is outside G4VisExtent", logger );
  } 

  for( i=0; i < n; i++ ) {
    G4ThreeVector vr = (*outside)[i] - point;
    G4ThreeVector v = vr.unit();

    G4bool validNorm;
    G4ThreeVector norm;

    G4double dist = testVolume->DistanceToOut( point, v, true, &validNorm, &norm );
    G4double NormalDist ;

    NormalDist = testVolume->DistanceToOut( point );

    if (dist <= 0) {
      ReportError( nError, point, v, safeDistance, "TI: DistanceToOut(p,v) <= 0  Normal Dist = ", logger );
      continue;
    }
    if (dist == kInfinity) {
      ReportError( nError, point, v, safeDistance, "TI: DistanceToOut(p,v) == kInfinity", logger );
      continue;
    }
    if (dist < safeDistance-1E-10) {
      ReportError( nError, point, v, safeDistance, "TI: DistanceToOut(p,v) < DistanceToIn(p)", logger );
      continue;
    }

    if (validNorm) {
      if (norm.dot(v) < 0) {
	ReportError( nError, point, v, safeDistance, "TI: Outgoing normal incorrect", logger );
	continue;
      }
    }

    G4ThreeVector p = point + v*dist;

    EInside insideOrNot = testVolume->Inside(p);
    if (insideOrNot == kInside) {
      ReportError( nError, point, v, safeDistance, "TI: DistanceToOut(p,v) undershoots", logger );
      continue;
    }
    if (insideOrNot == kOutside) {
      ReportError( nError, point, v, safeDistance, "TI: DistanceToOut(p,v) overshoots", logger );
      continue;
    }

    G4ThreeVector norm1 = testVolume->SurfaceNormal( p );
    if (norm1.dot(v) < 0) {
      if (testVolume->DistanceToIn(p,v) != 0)
	ReportError( nError, p, v, safeDistance, "TI: SurfaceNormal is incorrect", logger );
    }
  }
}


//
// ReportError
//
// Report the specified error message, but only if it has not been reported a zillion
// times already.
//
void SBTrun::ReportError( G4int *nError, const G4ThreeVector p, 
			  const G4ThreeVector v, const G4double distance,
			  const G4String comment, G4std::ostream &logger )
{
  //
  // Have we encountered this error message before?
  //
  SBTrunErrorList *last, *errors = errorList;
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
    errors = new SBTrunErrorList;
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
  logger << " DistanceToIn = " << distance ; 
  logger << G4endl;

  logger << ++(*nError) << " " << p.x() << " " << p.y() << " " << p.z() 
	 << " " << v.x() << " " << v.y() << " " << v.z() << G4endl;
}


//
// ClearErrors
//
// Reset list of errors (and clear memory)
//
void SBTrun::ClearErrors()
{
  SBTrunErrorList *here, *next;

  here = errorList;
  while( here ) {
    next = here->next;
    delete here;
    here = next;
  }
  errorList = 0;
}


//
// CountErrors
//
// Count up all errors
//
G4int SBTrun::CountErrors() const
{
  SBTrunErrorList *here;
  G4int answer = 0;

  here = errorList;
  while( here ) {
    answer += here->nUsed;
    here = here->next;
  }

  return answer;
}


//
// GetLoggedPV
//
// Get the p and v vectors stored in a test3 log file
//
G4int SBTrun::GetLoggedPV( G4std::istream &logger, const G4int errorIndex,
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

G4String SBTrun::CurrentSolid = "" ;

G4String SBTrun::GetCurrentSolid(void) 
{
  return CurrentSolid; 
}

void SBTrun::SetCurrentSolid(G4String SolidName)
{
  CurrentSolid = SolidName;
}

//
// --------------------------------------
// SBTrunPointList stuff
//

SBTrunPointList::SBTrunPointList( G4int size )
{
  pointList = new G4ThreeVector[size];
  maxPoints = size;
  numPoints = 0;
}

SBTrunPointList::~SBTrunPointList()
{
  delete [] pointList;
}

void SBTrunPointList::AddPoint( G4ThreeVector newPoint )
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
    G4int	irand = (G4int)((G4UniformRand()*( (G4double)maxPoints )));
    pointList[irand] = newPoint;
  }
}



