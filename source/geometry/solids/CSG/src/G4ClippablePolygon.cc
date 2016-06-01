//
// G4ClippablePolygon.cc
//
// Based on code from G4VSolid (P. Kent, V. Grichine, J. Allison)
//

#include "G4ClippablePolygon.hh"

#include "G4VoxelLimits.hh"

//
// AddVertexInOrder
//
void G4ClippablePolygon::AddVertexInOrder( const G4ThreeVector vertex )
{
	vertices.append( vertex );
}


//
// ClearAllVertices
//
void G4ClippablePolygon::ClearAllVertices()
{
	vertices.clear();
}


//
// Clip
//
void G4ClippablePolygon::Clip( const G4VoxelLimits &voxelLimit )
{    
	//
	// Heh. Do we have anything to do?
	//
	if (!voxelLimit.IsLimited()) return;
	
	//
	// Loop over all axes
	//
	static EAxis axes[3] = { kXAxis, kYAxis, kZAxis };
	
	EAxis *axis = axes;
	do {
		if (voxelLimit.IsLimited(*axis)) {
			G4ThreeVectorList tempPolygon;
		
			//
			// Build a "simple" voxelLimit that includes only the min extent
			// and apply this to our vertices, producing result in tempPolygon
			//
			G4VoxelLimits simpleLimit1;
			simpleLimit1.AddLimit( *axis, voxelLimit.GetMinExtent(*axis), kInfinity );
			ClipToSimpleLimits( vertices, tempPolygon, simpleLimit1 );

			//
			// If nothing is left from the above clip, we might as well return now
			// (but with an empty vertices)
			//
			if (tempPolygon.entries() == 0) {
				vertices.clear();
				return;
			}
			
			//
			// Now do the same, but using a "simple" limit that includes only the max extent.
			// Apply this to out tempPolygon, producing result in vertices.
			//
			G4VoxelLimits simpleLimit2;
			simpleLimit2.AddLimit( *axis, -kInfinity, voxelLimit.GetMaxExtent(*axis) );
			ClipToSimpleLimits( tempPolygon, vertices, simpleLimit2 );
			
			//
			// If nothing is left, return now
			//
			if (vertices.entries() == 0) return;
		}
	} while( ++axis < axes + sizeof(axes)/sizeof(EAxis) );
}



//
// GetExtent
//
void G4ClippablePolygon::GetExtent( const EAxis axis, 
				    G4double &min, G4double &max )
{
	G4int noLeft = vertices.entries();
	
	G4int i;
	for( i=0; i<noLeft; i++ ) {
		G4double component = vertices(i).operator()( axis );
		if (component < min )
			min = component;
		else if (component > max )
			max = component;
	}
}


// pVoxelLimits must be only limited along one axis, and either the maximum
// along the axis must be +kInfinity, or the minimum -kInfinity
void G4ClippablePolygon::ClipToSimpleLimits( G4ThreeVectorList& pPolygon,
				             G4ThreeVectorList& outputPolygon,
				             const G4VoxelLimits& pVoxelLimit   )
{
    G4int i;
    G4int noVertices=pPolygon.entries();
    G4ThreeVector vEnd,vStart;
    
    outputPolygon.clear();
    
    for (i=0;i<noVertices;i++)
	{
	    vStart=pPolygon(i);
	    if (i==noVertices-1)
		{
		    vEnd=pPolygon(0);
		}
	    else
		{
		    vEnd=pPolygon(i+1);
		}

	    if (pVoxelLimit.Inside(vStart))
		{
		    if (pVoxelLimit.Inside(vEnd))
			{
// vStart and vEnd inside -> output end point
			    outputPolygon.insert(vEnd);
			}
		    else
			{
// vStart inside, vEnd outside -> output crossing point
			    pVoxelLimit.ClipToLimits(vStart,vEnd);
			    outputPolygon.insert(vEnd);
			}
		    
		}
	    else
		{
		    if (pVoxelLimit.Inside(vEnd))
			{
// vStart outside, vEnd inside -> output inside section
			    pVoxelLimit.ClipToLimits(vStart,vEnd);
			    outputPolygon.insert(vStart);
			    outputPolygon.insert(vEnd);
			}
		    else
// Both point outside -> no output
			{
			}

		}
	}

}
