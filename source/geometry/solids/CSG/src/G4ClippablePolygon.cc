//
// G4ClippablePolygon.cc
//
// Based on code from G4VSolid (P. Kent, V. Grichine, J. Allison)
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
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
	if (!voxelLimit.IsLimited()) return;
	
	ClipAlongOneAxis( voxelLimit, kXAxis );
	ClipAlongOneAxis( voxelLimit, kYAxis );
	ClipAlongOneAxis( voxelLimit, kZAxis );
}


//
// PartialClip
//
// Clip, while ignoring the indicated axis
//
void G4ClippablePolygon::PartialClip( const G4VoxelLimits &voxelLimit, const EAxis IgnoreMe )
{
	if (!voxelLimit.IsLimited()) return;
	
	if (IgnoreMe != kXAxis) ClipAlongOneAxis( voxelLimit, kXAxis );
	if (IgnoreMe != kYAxis) ClipAlongOneAxis( voxelLimit, kYAxis );
	if (IgnoreMe != kZAxis) ClipAlongOneAxis( voxelLimit, kZAxis );
}


//
// GetExtent
//
G4bool G4ClippablePolygon::GetExtent( const EAxis axis, 
				      G4double &min, G4double &max )
{
	//
	// Okay, how many entries do we have?
	//
	G4int noLeft = vertices.entries();
	
	//
	// Return false if nothing is left
	//
	if (noLeft == 0) return false;
	
	//
	// Initialize min and max to our first vertex
	//
	min = max = vertices(0).operator()( axis );
	
	//
	// Compare to the rest
	//
	G4int i;
	for( i=1; i<noLeft; i++ ) {
		G4double component = vertices(i).operator()( axis );
		if (component < min )
			min = component;
		else if (component > max )
			max = component;
	}
	
	return true;
}



//
// Clip along just one axis, as specified in voxelLimit
//
void G4ClippablePolygon::ClipAlongOneAxis( const G4VoxelLimits &voxelLimit, const EAxis axis )
{    
	if (!voxelLimit.IsLimited(axis)) return;
	
	G4ThreeVectorList tempPolygon;

	//
	// Build a "simple" voxelLimit that includes only the min extent
	// and apply this to our vertices, producing result in tempPolygon
	//
	G4VoxelLimits simpleLimit1;
	simpleLimit1.AddLimit( axis, voxelLimit.GetMinExtent(axis), kInfinity );
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
	simpleLimit2.AddLimit( axis, -kInfinity, voxelLimit.GetMaxExtent(axis) );
	ClipToSimpleLimits( tempPolygon, vertices, simpleLimit2 );

	//
	// If nothing is left, return now
	//
	if (vertices.entries() == 0) return;
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
