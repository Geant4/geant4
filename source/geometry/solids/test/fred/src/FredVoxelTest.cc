//
// FredVoxelTest.cc
//
// Implementation of voxel solid tests
//
#include "FredVoxelTest.hh"
#include "G4VSolid.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Color.hh"
#include "G4Polyline.hh"
#include "G4VisAttributes.hh"

//
// Creator
//
FredVoxelTest::FredVoxelTest( )
{
	inverseTransform = transform.Inverse();
	
	test.solid = 0;
}

//
// Destructor
//
FredVoxelTest::~FredVoxelTest()
{;}


//
// SetExtent
//
void FredVoxelTest::SetExtent( const EAxis axis, const G4double min, const G4double max )
{
	voxelLimits.AddLimit( axis, min, max );
}


//
// SetOrigin
//
void FredVoxelTest::SetOrigin( const G4ThreeVector origin )
{
	transform.SetNetTranslation( origin );
	inverseTransform = transform.Inverse();
}


//
// Rotate
//	
void FredVoxelTest::Rotate( const EAxis axis, const G4double value )
{
	switch(axis) {
		case kXAxis: rotation = rotation.rotateX(value); break;
		case kYAxis: rotation = rotation.rotateY(value); break;
		case kZAxis: rotation = rotation.rotateZ(value); break;
	}
	transform.SetNetRotation( rotation );
	inverseTransform = transform.Inverse();
}


//
// Reset rotation
//
void FredVoxelTest::ResetRotation()
{
	rotation = G4RotationMatrix();
	transform.SetNetRotation( rotation );
	inverseTransform = transform.Inverse();
}


//
// Test
//
void FredVoxelTest::Test( const EAxis axis, const G4VSolid *solid )
{
	test.solid = solid;
	test.axis  = axis;
	
	test.result = solid->CalculateExtent( axis, voxelLimits, transform, test.min, test.max );
	if (test.result) {
		G4cout << "Voxel intersects the solid" << G4endl;
	}
	else {
		G4cout << "Voxel does not intersect the solid" << G4endl;
	}
}


//
// Draw
//
void FredVoxelTest::Draw( )
{
	static EAxis axes[3] = { kXAxis, kYAxis, kZAxis };
	static G4ThreeVector cartAxes[3] = { G4ThreeVector(1,0,0), 
				             G4ThreeVector(0,1,0),
				             G4ThreeVector(0,0,1) };
					     
        G4VVisManager *visManager = G4VVisManager::GetConcreteInstance();
	if (!visManager) return;

	G4ThreeVector origin(0,0,0);
	
	inverseTransform.ApplyPointTransform( origin );
	PlotMarker( origin, visManager );
	
	G4int i;
	for( i=0; i<3; i++ ) {
		if (!voxelLimits.IsLimited(axes[i])) continue;
		
		G4ThreeVector max(0,0,0);
		max = max + cartAxes[i]*voxelLimits.GetMaxExtent(axes[i]);
		
		inverseTransform.ApplyPointTransform( max );
		PlotMarker( max, visManager );
		PlotLine( origin, max, visManager );
		
		G4ThreeVector min(0,0,0);
		min = min + cartAxes[i]*voxelLimits.GetMinExtent(axes[i]);
		
		inverseTransform.ApplyPointTransform( min );
		PlotMarker( min, visManager );
		PlotLine( origin, min, visManager );
	}
	
	if (test.solid && test.result) {
		G4ThreeVector max, min, a, b, c, d;
	
		if (test.max < voxelLimits.GetMaxExtent(test.axis)) {
			max = cartAxes[test.axis]*test.max;
			a = inverseTransform.TransformPoint( max );
			PlotMarker( a, visManager, true );
			PlotLine( origin, a, visManager, true );
		}
	
		if (test.min > voxelLimits.GetMinExtent(test.axis)) {
			min = cartAxes[test.axis]*test.min;
			a = inverseTransform.TransformPoint( min );
			PlotMarker( a, visManager, true );
			PlotLine( origin, a, visManager, true );
		}
		
		for( i=0; i<3; i++ ) {
			if (axes[i]==test.axis) continue;
			
			if (voxelLimits.IsLimited(axes[i])) {
				G4ThreeVector aMax, bMax;
					      
				G4int j = (i + 1)%3;
				if (axes[j]==test.axis) j = (j + 1)%3;
				
				if (voxelLimits.IsLimited(axes[j])) {
					c = voxelLimits.GetMaxExtent(axes[j])*cartAxes[j];
					d = voxelLimits.GetMinExtent(axes[j])*cartAxes[j];
				}
				else {
					c = +4*m*cartAxes[j];
					d = -4*m*cartAxes[j];
				}
				
				if (test.max < voxelLimits.GetMaxExtent(test.axis)) {
					aMax = max + voxelLimits.GetMaxExtent(axes[i])*cartAxes[i];
					bMax = max + voxelLimits.GetMinExtent(axes[i])*cartAxes[i];

					a = inverseTransform.TransformPoint( aMax );
					b = inverseTransform.TransformPoint( bMax );
					PlotLine( a, b, visManager, true );

					a = inverseTransform.TransformPoint( aMax + c );
					b = inverseTransform.TransformPoint( aMax + d );
					PlotLine( a, b, visManager, true );
					a = inverseTransform.TransformPoint( bMax + c );
					b = inverseTransform.TransformPoint( bMax + d );
					PlotLine( a, b, visManager, true );
				}

				if (test.min > voxelLimits.GetMinExtent(test.axis)) {
					aMax = min + voxelLimits.GetMaxExtent(axes[i])*cartAxes[i];
					bMax = min + voxelLimits.GetMinExtent(axes[i])*cartAxes[i];

					a = inverseTransform.TransformPoint( aMax );
					b = inverseTransform.TransformPoint( bMax );
					PlotLine( a, b, visManager, true );

					a = inverseTransform.TransformPoint( aMax + c );
					b = inverseTransform.TransformPoint( aMax + d );
					PlotLine( a, b, visManager, true );
					a = inverseTransform.TransformPoint( bMax + c );
					b = inverseTransform.TransformPoint( bMax + d );
					PlotLine( a, b, visManager, true );
				}
			}
			else {
				G4ThreeVector aMax = max + 4*cartAxes[i],
					      bMax = max - 4*cartAxes[i];
					      
				if (test.max < voxelLimits.GetMaxExtent(test.axis)) {
					a = inverseTransform.TransformPoint( max + 4*m*cartAxes[i] );
					b = inverseTransform.TransformPoint( max - 4*m*cartAxes[i] );
					PlotLine( a, b, visManager, true );
				}
					      
				if (test.min > voxelLimits.GetMinExtent(test.axis)) {
					a = inverseTransform.TransformPoint( min + 4*m*cartAxes[i] );
					b = inverseTransform.TransformPoint( min - 4*m*cartAxes[i] );
					PlotLine( a, b, visManager, true );
				}
			}
		}
	}
}


//
// PlotMarker
//
void FredVoxelTest::PlotMarker( G4ThreeVector point, 
			        G4VVisManager *visManager, G4bool isAtest ) 
{
	G4Circle circle( point );
	circle.SetWorldSize( 5*cm );
	G4Color color( isAtest ? 1.0 : 0.25, isAtest ? 1.0 : 0.25, 1.0 );
	G4VisAttributes attribs( color );
	circle.SetVisAttributes( attribs );

	circle.SetFillStyle( G4Circle::filled );
	visManager->Draw( circle );
}


//
// PlotLine
//
void FredVoxelTest::PlotLine( G4ThreeVector start, G4ThreeVector end, 
			      G4VVisManager *visManager, G4bool isAtest ) 
{
	G4Polyline line;
	
	line.append( start );
	line.append( end );
	
	G4Color color( isAtest ? 1.0 : 0.25, isAtest ? 1.0 : 0.25, 1.0 );
	G4VisAttributes attribs( color );
	line.SetVisAttributes( attribs );
	
	visManager->Draw( line );
}
	


	
		
