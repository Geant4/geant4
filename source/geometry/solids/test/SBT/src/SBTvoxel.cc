//
// SBTvoxel.cc
//
// Implementation of a batch based voxel test
//

#include "globals.hh"
#include "Randomize.hh"

#include "SBTvoxel.hh"

#include "SBTVisManager.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

#include "G4VSolid.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include <time.h>
#include "g4std/iomanip"
#include "g4std/strstream"

//
// Constructor
//
SBTvoxel::SBTvoxel()
{
	SetDefaults();
}


//
// Destructor
//
SBTvoxel::~SBTvoxel() {;}


//
// SetDefaults
//
// Set default values for test parameters
//
void SBTvoxel::SetDefaults()
{
	target = G4ThreeVector( 0, 0, 0 );
	widths = G4ThreeVector( 1*m, 1*m, 1*m );
	
	maxVoxels = 100;
	maxErrors = 20;
}



//
// Debug
//
// Invoke the CalculateExtent method of the target volume.
// This can be particularly useful for debugging.
//
void SBTvoxel::Debug( const G4VSolid *testVolume, const EAxis axis,
		      const G4VoxelLimits &voxel, const G4AffineTransform &transform,
		      const G4ThreeVector *point ) const
{
	G4double min, max;
	
	if (testVolume->CalculateExtent( axis, voxel, transform, min, max ))
		G4cout << "Solid is intersected with: " << min << " " << max << G4endl;
	else
		G4cout << "Voxel misses solid" << G4endl;

	if (point) testVolume->Inside( *point );
}


//
// Draw
//
// Draw a test shape and a voxel with arbitrary limits, using
// an arbitrary tranformation, and with any number of optional
// markers
//
void SBTvoxel::Draw( const G4VSolid *testVolume, 
		     const G4VoxelLimits &voxel, const G4AffineTransform &transform,
		     const G4ThreeVector *points, const G4int numPoints,
		     const EAxis	axis, const G4double limits[2],
	 	     SBTVisManager *visManager ) const
{
	//
	// Get inverse transform
	//
	G4AffineTransform inverseTransform = transform.Inverse();

	//
	// Prepare visualization
	//
        visManager->ClearView();

	//
	// Draw voxel as a box.
	//
        G4VisAttributes blueStuff( G4Color(0,0,1) );
        G4VisAttributes yuckStuff( G4Color(0,1,0) );

	static const EAxis axes[3] = { kXAxis, kYAxis, kZAxis };
	static const G4ThreeVector axisVectors[3] = {G4ThreeVector(1,0,0),
						     G4ThreeVector(0,1,0),
						     G4ThreeVector(0,0,1) };

	G4bool drawLimits = limits[0] <= limits[1];

	G4ThreeVector dmin[3], dmax[3], lvec[2];
	G4int i;
	for( i=0; i<3; i++ ) {
		G4double min, max;
		if (voxel.IsLimited(axes[i])) {
			min = voxel.GetMinExtent(axes[i]);
			max = voxel.GetMaxExtent(axes[i]);
		}
		else {
			min = -4.0*m;
			max = +4.0*m;
		}

		dmin[i] = min*axisVectors[i];
		dmax[i] = max*axisVectors[i];
		
		if (drawLimits && axis==axes[i]) {
			lvec[0] = limits[0]*axisVectors[i];
			lvec[1] = limits[1]*axisVectors[i];
		}
	}
	
	for( i=0; i<3 ;i++ ) {
		G4int bitmask;
		for( bitmask=0; bitmask < 8; bitmask++ ) {
			if (bitmask&(1<<i)) continue;
			
			G4ThreeVector a = dmin[i], b = dmax[i];
			G4ThreeVector c = a, d = b;
			G4int j;
			for( j=0; j<3; j++ ) {
				if (i==j) continue;
				G4ThreeVector *use = bitmask&(1<<j) ? dmin+j : dmax+j;
				a += *use;
				b += *use;
				if (drawLimits) {
					if (axes[j]==axis) use = lvec + (bitmask&(1<<j) ? 0 : 1);
					c += *use;
					d += *use;
				}
			}
			G4Polyline polyline;
			polyline.SetVisAttributes( blueStuff );
			inverseTransform.ApplyPointTransform( a );
			inverseTransform.ApplyPointTransform( b );
			polyline.append( a );
			polyline.append( b );
			visManager->Draw( polyline );
			
			if (drawLimits && axes[i]!=axis) {
				G4Polyline polyline;
				polyline.SetVisAttributes( yuckStuff );
				inverseTransform.ApplyPointTransform( c );
				inverseTransform.ApplyPointTransform( d );
				polyline.append( c );
				polyline.append( d );
				visManager->Draw( polyline );
			}
		}
	}

	
	//
	// Draw points
	//
        G4VisAttributes whiteStuff( G4Color(1,1,1) );

	const G4ThreeVector *thisPoint;
	for (thisPoint = points; thisPoint < points+numPoints; thisPoint++ ) {
		G4Circle circle(*thisPoint);
		circle.SetWorldSize( 5*cm );
		circle.SetVisAttributes( whiteStuff );
		visManager->Draw( circle );
	}
		
	
	//
	// This draws the target solid
	//
        G4VisAttributes redStuff( G4Color(1,0,0) );
	visManager->Draw( *testVolume, redStuff );

        visManager->Show();;
}


//
// RunTest
//
// Perform a test on the specified solid
//
void SBTvoxel::RunTest( const G4VSolid *testVolume, G4std::ostream &logger )
{
	//
	// Output test parameters
	//
	time_t now = time(0);
	G4String dateTime( ctime(&now) );		// AFAIK, this is standard c++
	
	logger << "% SBT voxel logged output " << dateTime << G4endl;
	logger << "% target =    " << target << G4endl;
	logger << "% widths =    " << widths << G4endl;
	logger << "% maxVoxels = " << maxVoxels << G4endl;
	logger << "% maxErrors = " << maxErrors << G4endl;

	G4int nVoxel = 0,
	      nError = 0;
	      
	//
	// Generate a list of 1000 random points inside the solid
	//
	G4ThreeVector inside[1000];
	G4int numInside;
	
	GetInsidePoints( testVolume, inside, &numInside, 1000, 100000 );
	
	G4RotationMatrix randomRotate1, randomRotate2;
	      
	for(;;) {
		//
		// Generate a random voxel limit.
		// G4VoxelLimits has no "reset" method, so we need
		// to create a brand spanking new one each iteration
		//
		G4VoxelLimits *voxel = NewRandomVoxel( widths );

		//
		// Move it to random positions
		//
		G4int j;
		for( j=0; j < 20; j++ ) {
			G4ThreeVector offset = j > 0 ? GetRandomPoint() : G4ThreeVector(0,0,0);
		
			G4AffineTransform transform( offset );
			
			//
			// Test aligned
			//
			if (TestOneVoxel( testVolume, *voxel, transform, inside, numInside, logger )) {
				if (++nError >= maxErrors) break;
			}

			//
			// Generate a couple random orientations
			//
			randomRotate1.rotateZ( G4UniformRand() );
			transform.SetNetRotation( randomRotate1 );
			
			if (TestOneVoxel( testVolume, *voxel, transform, inside, numInside, logger )) {
				if (++nError >= maxErrors) break;
			}
			
			randomRotate2.rotateX( G4UniformRand() );
			transform.SetNetRotation( randomRotate2 );
			
			if (TestOneVoxel( testVolume, *voxel, transform, inside, numInside, logger )) {
				if (++nError >= maxErrors) break;
			}
			
			randomRotate2.rotateY( G4UniformRand() );
			transform.SetNetRotation( randomRotate2 );
			
			if (TestOneVoxel( testVolume, *voxel, transform, inside, numInside, logger )) {
				if (++nError >= maxErrors) break;
			}
			
			randomRotate2.rotateZ( G4UniformRand() );
			transform.SetNetRotation( randomRotate2 );
			
			if (TestOneVoxel( testVolume, *voxel, transform, inside, numInside, logger )) {
				if (++nError >= maxErrors) break;
			}
		}
		delete voxel;

		if (nError >= maxErrors) {
			logger << "% End of test (maximum number errors) ";
			break;
		}
		if (++nVoxel >= maxVoxels) {
			logger << "% End of test (maximum number voxels) ";
			break;
		}
	}
	
	now = time(0);
	G4String dateTime2( ctime(&now) );		
	logger << dateTime2 << G4endl;

	logger << "% Statistics: voxels=" << nVoxel << " errors=" << nError << G4endl;
	       
	logger << "%(End of file)" << G4endl;
}


//
// TestOneVoxel
//
G4bool SBTvoxel::TestOneVoxel( const G4VSolid *testVolume, 
			       const G4VoxelLimits &voxel,
			       const G4AffineTransform &transform,
			       const G4ThreeVector inside[], const G4int numInside,
			       G4std::ostream &logger ) const
{
	static const EAxis axes[3] = { kXAxis, kYAxis, kZAxis };
	G4int numError = 0;
	
	//
	// Get inverse transform
	//
	G4AffineTransform inverseTransform = transform.Inverse();
	
	//
	// Loop over the points, collecting min/max for each axis
	//
	G4double pointMins[3] = {+kInfinity, +kInfinity, +kInfinity}, 
		 pointMaxs[3] = {-kInfinity, -kInfinity, -kInfinity};
	G4int	 numPointInside = 0;
	
	G4int i = numInside;
	while( i-- > 0 ) {
		G4ThreeVector point = transform.TransformPoint( inside[i] );
		if (voxel.Inside(point)) {
			numPointInside++;
			
			G4int j;
			for (j=0; j<3; j++) {
				G4double pv = point(axes[j]);
				if (pv < pointMins[j]) pointMins[j] = pv;
				if (pv > pointMaxs[j]) pointMaxs[j] = pv;
			}
		}
	}

	//
	// Loop over axes
	//
	for( i=0; i<3; i++ ) {
		G4double min, max;
		
		//
		// Query the solid
		//
		if (testVolume->CalculateExtent( axes[i], voxel, transform, min, max )) {
			//
			// Compare min/max to the list of inside points
			//
			if (min > pointMins[i] || max < pointMaxs[i]) {
				numError++;
				logger << "ERROR: Voxel limits are incorrect, axis " << i << G4endl;
				logger << "  reported limits were: " << min << " " << max << G4endl;
				logger << "  but points were found at: " << pointMins[i] << " " << pointMaxs[i] << G4endl;
			}
			
			//
			// Min or max in reversee order?
			//
			if (min >= max) {
				numError++;
				logger << "ERROR: Voxel limits max <= min, axis " << i << G4endl;
				logger << "  reported limits were: " << min << " " << max << G4endl;
			}
			
			//
			// Min or max outside limits?
			//
			// We give the solid an extra "kCarTolerance" space, since
			// some solids like to add this value to their return values
			//
			if ( voxel.IsLimited(axes[i]) ) {
				if (min < voxel.GetMinExtent(axes[i])-1.1*kCarTolerance) {
					numError++;
					logger << "ERROR: Voxel min is below pre-existing limit, axis " << i << G4endl;
					logger << "  reported limits were: " << min << " " << max << G4endl;
				}
				if (max > voxel.GetMaxExtent(axes[i])+1.1*kCarTolerance) {
					numError++;
					logger << "ERROR: Voxel max is above pre-existing limit, axis " << i << G4endl;
					logger << "  reported limits were: " << min << " " << max << G4endl;
				}
			}
					
			
			if ( (!voxel.IsLimited(axes[i])) || max < voxel.GetMaxExtent(axes[i]) ) {
				//
				// Construct a set of points just outside the voxel limits
				// and make sure they are outside
				//
				G4ThreeVector testPoints[9];
				MakeVoxelTestPoints( voxel, axes[i], max, testPoints );
				
				G4ThreeVector *testPoint = testPoints;
				do {
					G4ThreeVector tp = inverseTransform.TransformPoint( *testPoint );
					if (testVolume->Inside( tp ) == kInside) {
						numError++;
						logger << "ERROR: Voxel MAX limit is too small, axis " << i << G4endl;
						logger << "  reported limits were: " << min << " " << max << G4endl;
						logger << "  test point " << *testPoint << " [" << tp << "] was inside" << G4endl;
						
						MakeVoxelTestPoints( voxel, axes[i], max+10.0*kCarTolerance, testPoints );
						tp = inverseTransform.TransformPoint( *testPoint );
						if (testVolume->Inside( tp ) == kInside)
							logger << "  and so was a more tolerant test point" << G4endl;
						break;
					}
				} while( ++testPoint < testPoints+9 );
			}
			
			if ( (!voxel.IsLimited(axes[i])) || min > voxel.GetMinExtent(axes[i]) ) {
				G4ThreeVector testPoints[9];
				MakeVoxelTestPoints( voxel, axes[i], min, testPoints );
				
				G4ThreeVector *testPoint = testPoints;
				do {
					G4ThreeVector tp = inverseTransform.TransformPoint( *testPoint );
					if (testVolume->Inside( tp ) == kInside) {
						numError++;
						logger << "ERROR: Voxel MIN limit is too large, axis " << i << G4endl;
						logger << "  reported limits were: " << min << " " << max << G4endl;
						logger << "  test point " << *testPoint << " [" << tp << "] was inside" << G4endl;
						
						MakeVoxelTestPoints( voxel, axes[i], min-10.0*kCarTolerance, testPoints );
						tp = inverseTransform.TransformPoint( *testPoint );
						if (testVolume->Inside( tp ) == kInside)
							logger << "  and so was a more tolerant test point" << G4endl;
						break;
					}
				} while( ++testPoint < testPoints+9 );
			}
			
		}
		else {
			//
			// The voxel does not intersect the solid
			// Make sure there are no inside points inside
			//
			if (numPointInside) {
				numError++;
				logger << "ERROR: Voxel isn't empty, axis " << i 
				       << ",  number points inside: " << numPointInside << " out of " << numInside << G4endl;
				logger << "  they have limits of: " << pointMins[i] << " " << pointMaxs[i] << G4endl;
			}
		}
	}
	
	if (numError) {
		DumpVoxel( voxel, logger );
		DumpTransform( transform, logger );
		return true;
	}
	
	return false;
}


//
// DumpVoxel
//
void SBTvoxel::DumpVoxel( const G4VoxelLimits &voxel, G4std::ostream &logger ) const
{
	logger << "VOXEL =";
	
	static const EAxis axes[3] = { kXAxis, kYAxis, kZAxis };
	
	const EAxis *axis = axes;
	do {
		logger << " (";
		if (voxel.IsLimited(*axis))
			logger << voxel.GetMinExtent(*axis) << " "
			       << voxel.GetMaxExtent(*axis);
		else
			logger << "unlimited";
		logger << ")";
	} while( ++axis < axes+3 );
	
	logger << G4endl;
}	


//
// DumpTransform
//
void SBTvoxel::DumpTransform( const G4AffineTransform &transform, G4std::ostream &logger ) const
{
	G4RotationMatrix rotate = transform.NetRotation();
  
  	G4ThreeVector axis;
	G4double      amount;
	
	rotate.getAngleAxis( amount, axis );
  
	logger << "TRANLATE = " << transform.NetTranslation() << G4endl;
	logger << "ROTATE = " << axis << " " << amount << G4endl;
}

	

//
// GetInsidePoints
//
void SBTvoxel::GetInsidePoints( const G4VSolid *testVolume,
				G4ThreeVector inside[], G4int *numInside,
				const G4int numPoints,
				const G4int maxAttempts ) const
{
	*numInside = 0;
	
	G4int i;
	for( i=0; i<maxAttempts; i++ ) {
		G4ThreeVector point = GetRandomPoint();
		
		if (testVolume->Inside( point ) == kInside) {
			inside[*numInside] = point;
			if (++(*numInside) == numPoints) return;
		}
	}
}


//
// GetRandomPoint
//
// Return a random point in three dimensions using the current
// specs 
//
G4ThreeVector SBTvoxel::GetRandomPoint() const {
	G4double dx = widths.x()*GaussianRandom(10*m/widths.x()),
		 dy = widths.y()*GaussianRandom(10*m/widths.y()),
		 dz = widths.z()*GaussianRandom(10*m/widths.z());
		 
		 
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
G4double SBTvoxel::GaussianRandom(const G4double cutoff) const {
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
// GetRandomLimit
//
G4bool SBTvoxel::GetRandomLimit( G4double x[2], const G4double range ) const
{
	//
	// Generate a random number
	//
	G4double rand = G4UniformRand();
	
	//
	// Let's say that, well, 20% of the time we are
	// not limited in this dimension
	//
	if (rand < 0.2) return false;
	
	//
	// Otherwise, construct limits
	//
	x[0] = range*(rand - 0.6)/0.4;
	x[1] = x[0] + range*G4UniformRand();
	return true;
}


//
// GetRandomVoxel
//
// Make a random voxel
//
G4VoxelLimits *SBTvoxel::NewRandomVoxel( const G4ThreeVector &theWidths ) const
{
	G4double xlim[2];
	
	G4VoxelLimits *voxel = new G4VoxelLimits;
	
	if (GetRandomLimit(xlim,widths.x())) voxel->AddLimit( kXAxis, xlim[0], xlim[1] );
	if (GetRandomLimit(xlim,widths.y())) voxel->AddLimit( kYAxis, xlim[0], xlim[1] );
	if (GetRandomLimit(xlim,widths.z())) voxel->AddLimit( kZAxis, xlim[0], xlim[1] );

	return voxel;
}
	

//
// MakeVoxelTestPoints
//
// Make a set of nine vectors that lay in a grid inside the
// voxel at the specified location along the specified axis
//
void SBTvoxel::MakeVoxelTestPoints( const G4VoxelLimits &voxel,
				    const EAxis valueAxis, const G4double value,
				    G4ThreeVector testPoints[9] ) const
{
	G4ThreeVector grid[6], offset;

	static const EAxis axes[3] = { kXAxis, kYAxis, kZAxis };
	static const G4ThreeVector axisVectors[3] = {G4ThreeVector(1,0,0),
						     G4ThreeVector(0,1,0),
						     G4ThreeVector(0,0,1) };
	
	const EAxis *axis = axes;
	G4ThreeVector *nextGrid = grid;
	const G4ThreeVector *axisVector = axisVectors;
	do {
		if (*axis == valueAxis) {
			offset = value*(*axisVector);
		}
		else if (voxel.IsLimited(*axis)) {
			G4double min = voxel.GetMinExtent(*axis);
			G4double max = voxel.GetMaxExtent(*axis);
			
			if (min <= -kInfinity) 
				min = max - 10.0*m;
			else if (max >= kInfinity)
				max = min + 10.0*m;
			
			(*nextGrid++) = 0.5*(min+max)*(*axisVector);
			(*nextGrid++) =           min*(*axisVector);
			(*nextGrid++) =           max*(*axisVector);
		}
		else {
			nextGrid++;		// zero vector
			(*nextGrid++) = +10.0*m*(*axisVector);
			(*nextGrid++) = -10.0*m*(*axisVector);
		}
	} while( ++axisVector, ++axis < axes+3 );
	
	testPoints[0] = offset + grid[0] + grid[3];
	testPoints[1] = offset + grid[1] + grid[3];
	testPoints[2] = offset + grid[2] + grid[3];
	testPoints[3] = offset + grid[0] + grid[4];
	testPoints[4] = offset + grid[1] + grid[4];
	testPoints[5] = offset + grid[2] + grid[4];
	testPoints[6] = offset + grid[0] + grid[5];
	testPoints[7] = offset + grid[1] + grid[5];
	testPoints[8] = offset + grid[2] + grid[5];
}
	
	

	
