//
// G4PolyhedraSide.cc
//
// Implemenation of the face representing one segmented side of a Polyhedra
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4PolyhedraSide.hh"
#include "G4IntersectingCone.hh"
#include "G4ClippablePolygon.hh"
#include "G4AffineTransform.hh"
#include "G4SolidExtentList.hh"

//
// Constructor
//
// Values for r1,z1 and r2,z2 should be specified in clockwise
// order in (r,z).
//
G4PolyhedraSide::G4PolyhedraSide( const G4PolyhedraSideRZ *prevRZ,
				  const G4PolyhedraSideRZ *tail,
				  const G4PolyhedraSideRZ *head,
				  const G4PolyhedraSideRZ *nextRZ,
				  const G4int theNumSide, 
				  const G4double thePhiStart, 
				  const G4double thePhiTotal, 
				  const G4bool thePhiIsOpen,
				  const G4bool isAllBehind )
{
	//
	// Record values
	//
	r[0] = tail->r; z[0] = tail->z;
	r[1] = head->r; z[1] = head->z;
	
	G4double phiTotal;
	
	phiIsOpen = thePhiIsOpen;
	if (phiIsOpen) {
		phiTotal = thePhiTotal;
		startPhi = thePhiStart;

		//
		// Set phi values to our conventions
		//
		while (startPhi < 0.0) startPhi += 2.0*M_PI;
	}
	else {
		phiTotal = 2*M_PI;
		startPhi = 0;
	}
	
	allBehind = isAllBehind;
		
	//
	// Make our intersecting cone
	//
	cone = new G4IntersectingCone( r, z );
	
	//
	// Construct side plane vector set
	//
	numSide = theNumSide;
	deltaPhi = phiTotal/theNumSide;
	endPhi = startPhi+phiTotal;
	
	vecs = new G4PolyhedraSideVec[numSide];
	
	edges = new G4PolyhedraSideEdge[phiIsOpen ? numSide+1 : numSide];
	
	//
	// ...this is where we start
	//
	G4double phi = startPhi;
	G4ThreeVector a1( r[0]*cos(phi), r[0]*sin(phi), z[0] ),
		      b1( r[1]*cos(phi), r[1]*sin(phi), z[1] ),
		      c1( prevRZ->r*cos(phi), prevRZ->r*sin(phi), prevRZ->z ),
		      d1( nextRZ->r*cos(phi), nextRZ->r*sin(phi), nextRZ->z ),
		      a2, b2, c2, d2;
	G4PolyhedraSideEdge *edge = edges;
		      
	G4PolyhedraSideVec *vec = vecs;
	do {
		//
		// ...this is where we are going
		//
		phi += deltaPhi;
		a2 = G4ThreeVector( r[0]*cos(phi), r[0]*sin(phi), z[0] );
		b2 = G4ThreeVector( r[1]*cos(phi), r[1]*sin(phi), z[1] );
		c2 = G4ThreeVector( prevRZ->r*cos(phi), prevRZ->r*sin(phi), prevRZ->z );
		d2 = G4ThreeVector( nextRZ->r*cos(phi), nextRZ->r*sin(phi), nextRZ->z );
		
		G4ThreeVector tt;	
		
		//
		// ...build some relevant vectors.
		//    the point is to sacrifice a little memory with precalcs 
		//    to gain speed
		//
		vec->center = 0.25*( a1 + a2 + b1 + b2 );
		
		tt = b2 + b1 - a2 - a1;
		vec->surfRZ = tt.unit();
		if (vec==vecs) lenRZ = 0.25*tt.mag();
		
		tt = b2 - b1 + a2 - a1;
		vec->surfPhi = tt.unit();
		if (vec==vecs) {
			lenPhi[0] = 0.25*tt.mag();
			tt = b2 - b1;
			lenPhi[1] = (0.5*tt.mag()-lenPhi[0])/lenRZ;
		}
		
		tt = vec->surfPhi.cross(vec->surfRZ);
		vec->normal = tt.unit();
		
		//
		// ...edge normals are the average of the normals of
		//    the two faces they connect.
		//
		// ...edge normals are necessary if we are to accurately
		//    decide if a point is "inside" a face. For non-convex
		//    shapes, it is absolutely necessary to know information
		//    on adjacent faces to accurate determine this.
		//
		// ...we don't need them for the phi edges, since that
		//    information is taken care of internally. The r/z edges,
		//    however, depend on the adjacent G4PolyhedraSide.
		//
		G4ThreeVector a12, adj;
		
		a12 = a2-a1;

		adj = 0.5*(c1+c2-a1-a2);
		adj = adj.cross(a12);	
		adj = adj.unit() + vec->normal;	     
		vec->edgeNorm[0] = adj.unit();
		
		a12 = b1-b2;
		adj = 0.5*(d1+d2-b1-b2);
		adj = adj.cross(a12);	
		adj = adj.unit() + vec->normal;	     
		vec->edgeNorm[1] = adj.unit();
		
		//
		// ...the corners are crucial. It is important that
		//    they are calculated consistently for adjacent
		//    G4PolyhedraSides, to avoid gaps caused by roundoff.
		//
		vec->edges[0] = edge;
		edge->corner[0] = a1;
		edge->corner[1] = b1;
		edge++;
		vec->edges[1] = edge;

		a1 = a2;
		b1 = b2;
		c1 = c2;
		d1 = d2;
	} while( ++vec < vecs+numSide );
	
	//
	// Clean up hanging edge
	//
	if (phiIsOpen) {
		edge->corner[0] = a2;
		edge->corner[1] = b2;
	}
	else {
		vecs[numSide-1].edges[1] = edges;
	}
	
	//
	// Go back and fill in remaining fields in edges
	//
	vec = vecs;
	G4PolyhedraSideVec *prev = vecs+numSide-1;
	do {
		edge = vec->edges[0];		// The edge between prev and vec
		
		//
		// Okay: edge normal is average of normals of adjacent faces
		//
		G4ThreeVector eNorm = vec->normal + prev->normal;
		edge->normal = eNorm.unit();	
		
		//
		// Vertex normal is average of norms of attached edges 
		//
		eNorm = edge->normal + vec->edgeNorm[0] + prev->edgeNorm[0];
		edge->cornNorm[0] = eNorm.unit();
	
		eNorm = edge->normal + vec->edgeNorm[1] + prev->edgeNorm[1];
		edge->cornNorm[1] = eNorm.unit();
	} while( prev=vec, ++vec < vecs + numSide );
	
	if (phiIsOpen) {
		G4double rFact = cos(0.5*deltaPhi);
		//
		// If phi is open, we need to patch up the first and last edges
		//
		G4double phi1 = startPhi - 0.5*M_PI;
		G4ThreeVector phiNorm( cos(phi1), sin(phi1), 0 );
		
		vec = vecs;
		
		//
		// Edge normal is average of vec->normal and the normal
		// of the face closing the polyhedra in phi
		//
		G4ThreeVector eNorm = vec->normal + phiNorm;
		vec->edges[0]->normal = eNorm.unit();
		
		//
		// We need the edge normals (like above) of the adjacent
		// G4PolyhedraSides.
		//
		G4double dr = r[0]-prevRZ->r, dz = z[0]-prevRZ->z;
		phi1 = startPhi + 0.5*deltaPhi;
		eNorm = G4ThreeVector( dz*rFact*cos(phi1), dz*rFact*sin(phi1), -dr );
		
		//
		// Average three line normals for the vertex normal
		//
		eNorm = eNorm.unit() + vec->edges[0]->normal + vec->edgeNorm[0];
		vec->edges[0]->cornNorm[0] = eNorm.unit();
		
		//
		// Repeat for edgeNorm[1]
		//
		dr = nextRZ->r-r[1], dz = nextRZ->z-z[1];
		eNorm = G4ThreeVector( dz*rFact*cos(phi1), dz*rFact*sin(phi1), -dr );
		eNorm = eNorm.unit() + vec->edges[0]->normal + vec->edgeNorm[1];
		vec->edges[0]->cornNorm[1] = eNorm.unit();
		
		//
		// That was bad...
		//
		// But, now repeat for ending phi (edge[1])
		//
		phi1 = startPhi + phiTotal + 0.5*M_PI;
		phiNorm = G4ThreeVector( cos(phi1), sin(phi1), 0 );
		
		vec = vecs + numSide - 1;

		eNorm = vec->normal + phiNorm;
		vec->edges[1]->normal = eNorm.unit();

		dr = r[0]-prevRZ->r, dz = z[0]-prevRZ->z;
		phi1 = startPhi + phiTotal - 0.5*deltaPhi;
		eNorm = G4ThreeVector( dz*rFact*cos(phi1), dz*rFact*sin(phi1), -dr );

		eNorm = eNorm.unit() + vec->edges[1]->normal + vec->edgeNorm[0];
		vec->edges[1]->cornNorm[0] = eNorm.unit();

		dr = nextRZ->r-r[1], dz = nextRZ->z-z[1];
		eNorm = G4ThreeVector( dz*rFact*cos(phi1), dz*rFact*sin(phi1), -dr );
		eNorm = eNorm.unit() + vec->edges[1]->normal + vec->edgeNorm[1];
		vec->edges[1]->cornNorm[1] = eNorm.unit();
		
		//
		// Phew! I need a beer!
		//
	}
	
	//
	// edgeNorm is the factor one multiplies the distance along vector phi
	// on the surface of one of our sides in order to calculate the distance
	// from the edge. (see routine DistanceAway)
	//
	edgeNorm = 1.0/sqrt( 1.0 + lenPhi[1]*lenPhi[1] );
}


//
// Destructor
//	
G4PolyhedraSide::~G4PolyhedraSide()
{
	delete cone;
	delete [] vecs;
	delete [] edges;
}


//
// Copy constructor
//
G4PolyhedraSide::G4PolyhedraSide( const G4PolyhedraSide &source )
{
	CopyStuff( source );
}


//
// Assignment operator
//
G4PolyhedraSide *G4PolyhedraSide::operator=( const G4PolyhedraSide &source )
{
	if (this == &source) return this;
	
	delete cone;
	delete [] vecs;
	delete [] edges;
	
	CopyStuff( source );

	return this;
}


//
// CopyStuff
//
void G4PolyhedraSide::CopyStuff( const G4PolyhedraSide &source )
{
	//
	// The simple stuff
	//
	numSide		= source.numSide;
	r[0]		= source.r[0];
	r[1]		= source.r[1];
	z[0]		= source.z[0];
	z[1]		= source.z[1];
	startPhi	= source.startPhi;
	deltaPhi	= source.deltaPhi;
	endPhi		= source.endPhi;
	phiIsOpen	= source.phiIsOpen;
	allBehind	= source.allBehind;
	
	lenRZ		= source.lenRZ;
	lenPhi[0]	= source.lenPhi[0];
	lenPhi[1]	= source.lenPhi[1];
	edgeNorm	= source.edgeNorm;
	
	cone = new G4IntersectingCone( *source.cone );

	//
	// Duplicate edges
	//
	G4int	numEdges = phiIsOpen ? numSide+1 : numSide;
	edges = new G4PolyhedraSideEdge[numEdges];
	
	G4PolyhedraSideEdge *edge = edges,
			    *sourceEdge = source.edges;
	do {
		*edge = *sourceEdge;
	} while( ++sourceEdge, ++edge < edges + numEdges);

	//
	// Duplicate vecs
	//
	vecs = new G4PolyhedraSideVec[numSide];
	
	G4PolyhedraSideVec *vec = vecs,
			   *sourceVec = source.vecs;
	do {
		*vec = *sourceVec;
		vec->edges[0] = edges + (sourceVec->edges[0] - source.edges);
		vec->edges[1] = edges + (sourceVec->edges[1] - source.edges);
	} while( ++sourceVec, ++vec < vecs + numSide );
}
	

//
// Intersect
//
// Decide if a line intersects the face.
//
// Arguments:
//	p		= (in) starting point of line segment
//	v		= (in) direction of line segment (assumed a unit vector)
//	A, B		= (in) 2d transform variables (see note top of file)
//	normSign	= (in) desired sign for dot product with normal (see below)
//	surfTolerance	= (in) minimum distance from the surface
//	vecs		= (in) Vector set array
//	distance	= (out) distance to surface furfilling all requirements
//	distFromSurface = (out) distance from the surface
//	thisNormal	= (out) normal vector of the intersecting surface
//
// Return value:
//	true if an intersection is found. Otherwise, output parameters are undefined.
//
// Notes:
//    * normSign: if we are "inside" the shape and only want to find out how far
//      to leave the shape, we only want to consider intersections with surfaces in
//      which the trajectory is leaving the shape. Since the normal vectors to the
//      surface always point outwards from the inside, this means we want the dot
//      product of the trajectory direction v and the normal of the side normals[i]
//      to be positive. Thus, we should specify normSign as +1.0. Otherwise, if
//      we are outside and want to go in, normSign should be set to -1.0.
//      Don't set normSign to zero, or you will get no intersections!
//
//    * surfTolerance: see notes on argument "surfTolerance" in routine "IntersectSidePlane".
//      ----HOWEVER---- We should *not* apply this surface tolerance if the starting
//      point is not within phi or z of the surface. Specifically, if the starting
//      point p angle in x/y places it on a separate side from the intersection or
//      if the starting point p is outside the z bounds of the segment, surfTolerance
//      must be ignored or we should *always* accept the intersection! 
//      This is simply because the sides do not have infinite extent.
//      
//
G4bool G4PolyhedraSide::Intersect( const G4ThreeVector &p, const G4ThreeVector &v,	
				   const G4bool outgoing, const G4double surfTolerance,
				   G4double &distance, G4double &distFromSurface,
				   G4ThreeVector &normal, G4bool &isAllBehind )
{
	G4double normSign = outgoing ? +1 : -1;
	
	//
	// ------------------TO BE IMPLEMENTED---------------------
	// Testing the intersection of individual phi faces is
	// pretty straight forward. The simple thing therefore is to
	// form a loop and check them all in sequence.
	//
	// But, I worry about one day someone making
	// a polygon with a thousands sides. A linear search
	// would not be ideal in such a case.
	//
	// So, it would be nice to be able to quickly decide
	// which face would be intersected. One can make a very
	// good guess by using the intersection with a cone.
	// However, this is only reliable in 99% of the cases.
	//
	// My solution: make a decent guess as to the one or
	// two potential faces might get intersected, and then
	// test them. If we have the wrong face, use the test
	// to make a better guess.
	//
	// Since we might have two guesses, form a queue of
	// potential intersecting faces. Keep an array of 
	// already tested faces to avoid doing one more than
	// once.
	//
	// Result: at worst, an iterative search. On average,
	// a little more than two tests would be required.
	//
	G4ThreeVector q = p + v;
	
	G4int face = 0;
	G4PolyhedraSideVec *vec = vecs;
	do {
		//
		// Correct normal?
		//
		G4double dotProd = normSign*v.dot(vec->normal);
		if (dotProd <= 0) continue;
	
		//
		// Is this face in front of the point along the trajectory?
		//
		G4ThreeVector delta = p - vec->center;
		distFromSurface = -normSign*delta.dot(vec->normal);
		
		if (distFromSurface < -surfTolerance) continue;
		
		//
		//                            phi
		//      c -------- d           ^
		//      |          |           |
		//      a -------- b           +---> r/z
		//
		//
		// Do we remain on this particular segment?
		//
		G4ThreeVector qc = q - vec->edges[1]->corner[0];
		G4ThreeVector qd = q - vec->edges[1]->corner[1];
		
		if (normSign*qc.cross(qd).dot(v) < 0) continue;
		
		G4ThreeVector qa = q - vec->edges[0]->corner[0];
		G4ThreeVector qb = q - vec->edges[0]->corner[1];
		
		if (normSign*qa.cross(qb).dot(v) > 0) continue;
		
		//
		// We found the one and only segment we might be intersecting.
		// Do we remain within r/z bounds?
		//
		
		if (normSign*qa.cross(qc).dot(v) < 0) return false;
		if (normSign*qb.cross(qd).dot(v) > 0) return false;
		
		//
		// We allow the face to be slightly behind the trajectory
		// (surface tolerance) only if the point p is within
		// the vicinity of the face
		//
		if (distFromSurface < 0) {
			G4ThreeVector ps = p - vec->center; 
			
			G4double rz = ps.dot(vec->surfRZ);
			if (fabs(rz) > lenRZ+surfTolerance) return false; 

			G4double pp = ps.dot(vec->surfPhi);
			if (fabs(pp) > lenPhi[0] + lenPhi[1]*rz + surfTolerance) return false;
		}
			

		//
		// Intersection found. Return answer.
		//
		distance = distFromSurface/dotProd;
		normal = vec->normal;
		isAllBehind = allBehind;
		return true;
	} while( ++vec, ++face < numSide );

	//
	// Oh well. Better luck next time.
	//
	return false;
}


G4double G4PolyhedraSide::Distance( const G4ThreeVector &p, const G4bool outgoing )
{
	G4double normSign = outgoing ? -1 : +1;
	
	//
	// Try the closest phi segment first
	//
	G4int iPhi = ClosestPhiSegment( p.phi() );
	
	G4ThreeVector pdotc = p - vecs[iPhi].center;
	G4double normDist = pdotc.dot(vecs[iPhi].normal);
	
	if (normSign*normDist > -0.5*kCarTolerance) {
		return DistanceAway( p, vecs[iPhi], &normDist );
	}

	//
	// Now we have an interesting problem... do we try to find the
	// closest facing side??
	//
	// Considered carefully, the answer is no. We know that if we
	// are asking for the distance out, we are supposed to be inside,
	// and vice versa.
	//
	
	return kInfinity;
}


//
// Inside
//
EInside G4PolyhedraSide::Inside( const G4ThreeVector &p, const G4double tolerance, 
				 G4double *bestDistance )
{
	//
	// Which phi segment is closest to this point?
	//
	G4int iPhi = ClosestPhiSegment( p.phi() );
	
	G4double norm;
	
	//
	// Get distance to this segment
	//
	*bestDistance = DistanceToOneSide( p, vecs[iPhi], &norm );
	
	//
	// Use distance along normal to decide return value
	//
	if ((fabs(norm) < tolerance) && (*bestDistance < 2.0*tolerance) )
		return kSurface;
	else if (norm < 0)
		return kInside;
	else	
		return kOutside;
}


//
// Normal
//
G4ThreeVector G4PolyhedraSide::Normal( const G4ThreeVector &p,  G4double *bestDistance )
{
	G4int iPhi = ClosestPhiSegment( p.phi() );
	return vecs[iPhi].normal;
}


//
// Extent
//
G4double G4PolyhedraSide::Extent( const G4ThreeVector axis )
{
	if (axis.perp2() < DBL_MIN) {
		//
		// Special case
		//
		return axis.z() < 0 ? -cone->ZLo() : cone->ZHi();
	}

	G4int iPhi, i1, i2;
	G4double best;
	G4ThreeVector *list[4];
	
	//
	// Which phi segment, if any, does the axis belong to
	//
	iPhi = PhiSegment( axis.phi() );
	
	if (iPhi < 0) {
		//
		// No phi segment? Check front edge of first side and
		// last edge of second side
		//
		i1 = 0; i2 = numSide-1;
	}
	else {
		//
		// Check all corners of matching phi side
		//
		i1 = iPhi; i2 = iPhi;
	}
	
	list[0] = vecs[i1].edges[0]->corner;
	list[1] = vecs[i1].edges[0]->corner+1;
	list[2] = vecs[i2].edges[1]->corner;
	list[3] = vecs[i2].edges[1]->corner+1;
				
	//
	// Who's biggest?
	//
	best = -kInfinity;
	G4ThreeVector **vec = list;
	do {
		G4double answer = (*vec)->dot(axis);
		if (answer > best) best = answer;
	} while( ++vec < list+4 );
	
	return best;
}


//
// CalculateExtent
//
// See notes in G4VCSGface
//
void G4PolyhedraSide::CalculateExtent( const EAxis axis, 
				       const G4VoxelLimits &voxelLimit,
				       const G4AffineTransform &transform,
				       G4SolidExtentList &extentList        )
{
	G4ClippablePolygon polygon;
	
	//
	// Loop over all sides
	//
	G4PolyhedraSideVec *vec = vecs;
	do {
		//
		// Fill our polygon with the four corners of
		// this side, after the specified transformation
		//
		polygon.ClearAllVertices();
		
		polygon.AddVertexInOrder( transform.TransformPoint( vec->edges[0]->corner[0] ) );
		polygon.AddVertexInOrder( transform.TransformPoint( vec->edges[0]->corner[1] ) );
		polygon.AddVertexInOrder( transform.TransformPoint( vec->edges[1]->corner[1] ) );
		polygon.AddVertexInOrder( transform.TransformPoint( vec->edges[1]->corner[0] ) );
		
		//
		// Get extent
		//	
		polygon.Clip( voxelLimit );
		
		//
		// Add it to the list
		//
		G4double min, max;
		
		if (polygon.GetExtent( axis, min, max )) {
			//
			// Get dot product of normal along target axis
			//
			G4ThreeVector rotatedNormal = transform.TransformAxis(vec->normal);
			G4double dotNormal = rotatedNormal(axis);
			
			extentList.AddSurface( min, max, dotNormal, dotNormal, kCarTolerance );
		}
	} while( ++vec < vecs+numSide );
	
	return;
}


//
// -------------------------------------------------------

//
// IntersectSidePlane
//
// Decide if a line correctly intersects one side plane of our segment.
// It is assumed that the correct side has been chosen, and thus only 
// the z bounds (of the entire segment) are checked.
//
// normSign - To be multiplied against normal:
//            = +1.0 normal is unchanged
//            = -1.0 normal is reversed (now points inward)
//
// Arguments:
//	p		- (in) Point
//	v		- (in) Direction
//	vec		- (in) Description record of the side plane
//	normSign	- (in) Sign (+/- 1) to apply to normal
//	surfTolerance	- (in) Surface tolerance (generally > 0, see below)
//	distance	- (out) Distance along v to intersection
//	distFromSurface - (out) Distance from surface normal
//
// Notes:
// 	surfTolerance	- Used to decide if a point is behind the surface,
//			  a point is allow to be -surfTolerance behind the
//			  surface (as measured along the normal), but *only*
//			  if the point is within the r/z bounds + surfTolerance
//			  of the segment.
//
G4bool G4PolyhedraSide::IntersectSidePlane( const G4ThreeVector &p, const G4ThreeVector &v,
					    const G4PolyhedraSideVec vec,
					    const G4double normSign, 
					    const G4double surfTolerance,
				 	    G4double &distance, G4double &distFromSurface )
{
	//
	// Correct normal? Here we have straight sides, and can safely ignore
	// intersections where the dot product with the normal is zero.
	//
	G4double dotProd = normSign*v.dot(vec.normal);
	
	if (dotProd <= 0) return false;
	
	//
	// Calculate distance to surface. If the side is too far
	// behind the point, we must reject it.
	//
	G4ThreeVector delta = p - vec.center;
	distFromSurface = -normSign*delta.dot(vec.normal);
		
	if (distFromSurface < -surfTolerance) return false;

	//
	// Calculate precise distance to intersection with the side
	// (along the trajectory, not normal to the surface)
	//
	distance = distFromSurface/dotProd;
	
	//
	// Do we fall off the r/z extent of the segment?
	//
	// Calculate this very, very carefully! Why?
	//         1. If a RZ end is at R=0, you can't miss!
	//         2. If you just fall off in RZ, the answer must
	//            be consistent with adjacent G4PolyhedraSide faces.
	// (2) implies that only variables used by other G4PolyhedraSide
	// faces may be used, which includes only: p, v, and the edge corners.
	// It also means that one side is a ">" or "<", which the other
	// must be ">=" or "<=". Fortunately, this isn't a new problem.
	// The solution below I borrowed from Joseph O'Rourke,
	// "Computational Geometry in C (Second Edition)"
	// See: http://cs.smith.edu/~orourke/
	//
	G4ThreeVector ic = p + distance*v - vec.center;
	G4double atRZ = vec.surfRZ.dot(ic);
	
	if (atRZ < 0) {
		if (r[0]==0) return true;		// Can't miss!
		
		if (atRZ < -lenRZ*1.2) return false;	// Forget it! Missed by a mile.
		
		G4ThreeVector q = p + v;		
		G4ThreeVector qa = q - vec.edges[0]->corner[0],
			      qb = q - vec.edges[1]->corner[0];
		G4ThreeVector qacb = qa.cross(qb);
		if (normSign*qacb.dot(v) < 0) return false;
		
		if (distFromSurface < 0) {
			if (atRZ < -lenRZ-surfTolerance) return false;
		}
	}
	else if (atRZ > 0) {
		if (r[1]==0) return true;		// Can't miss!
		
		if (atRZ > lenRZ*1.2) return false;	// Missed by a mile
		
		G4ThreeVector q = p + v;		
		G4ThreeVector qa = q - vec.edges[0]->corner[1],
			      qb = q - vec.edges[1]->corner[1];
		G4ThreeVector qacb = qa.cross(qb);
		if (normSign*qacb.dot(v) >= 0) return false;
		
		if (distFromSurface < 0) {
			if (atRZ > lenRZ+surfTolerance) return false;
		}
	}

	return true;
}


//
// LineHitsSegments
//
// Calculate which phi segments a line intersects in three dimensions.
// No check is made as to whether the intersections are within the z bounds of
// the segment.
//
G4int G4PolyhedraSide::LineHitsSegments( const G4ThreeVector &p, const G4ThreeVector &v,
				 	 G4int *i1, G4int *i2 )
{
	G4double s1, s2;
	//
	// First, decide if and where the line intersects the cone
	//
	G4int n = cone->LineHitsCone( p, v, &s1, &s2 );
	
	if (n==0) return 0;
	
	//
	// Try first intersection.
	//
	*i1 = PhiSegment( atan2( p.y() + s1*v.y(), p.x() + s1*v.x() ) );
	if (n==1) {
		return (*i1 < 0) ? 0 : 1;
	}
	
	//
	// Try second intersection
	//
	*i2 = PhiSegment( atan2( p.y() + s2*v.y(), p.x() + s2*v.x() ) );
	if (*i1 == *i2) return 0;
	
	if (*i1 < 0) {
		if (*i2 < 0) return 0;
		*i1 = *i2;
		return 1;
	}

	if (*i2 < 0) return 1;
	
	return 2;
}


//
// ClosestPhiSegment
//
// Decide which phi segment is closest in phi to the point.
// The result is the same as PhiSegment if there is no phi opening.
//
G4int G4PolyhedraSide::ClosestPhiSegment( const G4double phi0 )
{
	G4int iPhi = PhiSegment( phi0 );
	if (iPhi >= 0) return iPhi;
	
	//
	// Boogers! The points falls inside the phi segment.
	// Look for the closest point: the start, or  end
	//
	G4double phi = phi0;
	
	while( phi < startPhi ) phi += 2*M_PI;
	G4double d1 = phi-endPhi;

	while( phi > startPhi ) phi -= 2*M_PI;
	G4double d2 = startPhi-phi;
	
	return (d2 < d1) ? 0 : numSide-1;
}


//
// PhiSegment
//
// Decide which phi segment an angle belongs to, counting from zero.
// A value of -1 indicates that the phi value is outside the shape
// (only possible if phiTotal < 360 degrees).
//
G4int G4PolyhedraSide::PhiSegment( const G4double phi0 )
{
	//
	// How far are we from phiStart? Come up with a positive answer
	// that is less than 2*PI
	//
	G4double phi = phi0 - startPhi;
	while( phi < 0      ) phi += 2*M_PI;
	while( phi > 2*M_PI ) phi -= 2*M_PI;

	//
	// Divide
	//
	G4int answer = phi/deltaPhi;
	
	if (answer >= numSide) {
		if (phiIsOpen) {
			return -1;	// Looks like we missed
		}
		else {
			answer = numSide-1;	// Probably just roundoff
		}
	}
	
	return answer;
}


//
// DistanceToOneSide
//
// Arguments:
//	p	 - (in) Point to check
//	vec	 - (in) vector set of this side
//	normDist - (out) distance normal to the side or edge, as appropriate, signed
// Return value = total distance from the side
//
G4double G4PolyhedraSide::DistanceToOneSide( const G4ThreeVector &p,
					     const G4PolyhedraSideVec vec,
					     G4double *normDist )
{
	G4ThreeVector pc = p - vec.center;
	
	//
	// Get normal distance
	//
	*normDist = vec.normal.dot(pc);

	//
	// Add edge penalty
	//
	return DistanceAway( p, vec, normDist );
}


//
// DistanceAway
//
// Add distance from side edges, if necesssary, to total distance,
// and updates normDist appropriate depending on edge normals.
//
G4double G4PolyhedraSide::DistanceAway( const G4ThreeVector &p,
					const G4PolyhedraSideVec vec,
					G4double *normDist )
{
	G4double distOut2;
	G4ThreeVector pc = p - vec.center;
	G4double distFaceNorm = *normDist;
	
	//
	// Okay, are we inside bounds?
	//
	G4double pcDotRZ  = pc.dot(vec.surfRZ);
	G4double pcDotPhi = pc.dot(vec.surfPhi);
	
	//
	// Go through all permutations.
	//                                                   Phi
	//               |              |                     ^
	//           B   |      H       |   E                 |
	//        ------[1]------------[3]-----               |
	//               |XXXXXXXXXXXXXX|                     +----> RZ
	//           C   |XXXXXXXXXXXXXX|   F
	//               |XXXXXXXXXXXXXX|
	//        ------[0]------------[2]----
	//           A   |      G       |   D
	//               |              |
	//
	// It's real messy, but at least it's quick
	//
	
	if (pcDotRZ < -lenRZ) {
		G4double lenPhiZ = lenPhi[0] - lenRZ*lenPhi[1];
		G4double distOutZ = pcDotRZ+lenRZ;
		//
		// Below in RZ
		//
		if (pcDotPhi < -lenPhiZ) {
			//
			// ...and below in phi. Find distance to point (A)
			//
			G4double distOutPhi = pcDotPhi+lenPhiZ;
			distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
			G4ThreeVector pa = p - vec.edges[0]->corner[0];
			*normDist = pa.dot(vec.edges[0]->cornNorm[0]);
		}
		else if (pcDotPhi > lenPhiZ) {
			//
			// ...and above in phi. Find distance to point (B)
			//
			G4double distOutPhi = pcDotPhi-lenPhiZ;
			distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
			G4ThreeVector pb = p - vec.edges[1]->corner[0];
			*normDist = pb.dot(vec.edges[1]->cornNorm[0]);
		}
		else {
			//
			// ...and inside in phi. Find distance to line (C)
			//
			G4ThreeVector pa = p - vec.edges[0]->corner[0];
			distOut2 = distOutZ*distOutZ;
			*normDist = pa.dot(vec.edgeNorm[0]);
		}
	}
	else if (pcDotRZ > lenRZ) {
		G4double lenPhiZ = lenPhi[0] + lenRZ*lenPhi[1];
		G4double distOutZ = pcDotRZ-lenRZ;
		//
		// Above in RZ
		//
		if (pcDotPhi < -lenPhiZ) {
			//
			// ...and below in phi. Find distance to point (D)
			//
			G4double distOutPhi = pcDotPhi+lenPhiZ;
			distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
			G4ThreeVector pd = p - vec.edges[0]->corner[1];
			*normDist = pd.dot(vec.edges[0]->cornNorm[1]);
		}
		else if (pcDotPhi > lenPhiZ) {
			//
			// ...and above in phi. Find distance to point (E)
			//
			G4double distOutPhi = pcDotPhi-lenPhiZ;
			distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
			G4ThreeVector pe = p - vec.edges[1]->corner[1];
			*normDist = pe.dot(vec.edges[1]->cornNorm[1]);
		}
		else {
			//
			// ...and inside in phi. Find distance to line (F)
			//
			distOut2 = distOutZ*distOutZ;
			G4ThreeVector pd = p - vec.edges[0]->corner[1];
			*normDist = pd.dot(vec.edgeNorm[1]);
		}
	}
	else {
		G4double lenPhiZ = lenPhi[0] + pcDotRZ*lenPhi[1];
	 	//
		// We are inside RZ bounds
		// 
		if (pcDotPhi < -lenPhiZ) {
			//
			// ...and below in phi. Find distance to line (G)
			//
			G4double distOut = edgeNorm*(pcDotPhi+lenPhiZ);
			distOut2 = distOut*distOut;
			G4ThreeVector pd = p - vec.edges[0]->corner[1];
			*normDist = pd.dot(vec.edges[0]->normal);
		}
		else if (pcDotPhi > lenPhiZ) {
			//
			// ...and above in phi. Find distance to line (H)
			//
			G4double distOut = edgeNorm*(pcDotPhi-lenPhiZ);
			distOut2 = distOut*distOut;
			G4ThreeVector pe = p - vec.edges[1]->corner[1];
			*normDist = pe.dot(vec.edges[1]->normal);
		}
		else {
			//
			// Inside bounds! No penalty.
			//
			return fabs(distFaceNorm);
		}
	}
	return sqrt( distFaceNorm*distFaceNorm + distOut2 );
}
