// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PolyconeSide.cc,v 1.1 2000-04-07 11:02:07 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4PolyconeSide.cc
//
// Implementation of the face representing one conical side of a polycone
//
// --------------------------------------------------------------------

#include "G4PolyconeSide.hh"
#include "G4IntersectingCone.hh"
#include "G4ClippablePolygon.hh"
#include "G4AffineTransform.hh"
#include "meshdefs.hh"
#include "G4SolidExtentList.hh"

//
// Constructor
//
// Values for r1,z1 and r2,z2 should be specified in clockwise
// order in (r,z).
//
G4PolyconeSide::G4PolyconeSide( const G4PolyconeSideRZ *prevRZ,
                        	const G4PolyconeSideRZ *tail,
                        	const G4PolyconeSideRZ *head,
                        	const G4PolyconeSideRZ *nextRZ,
				const G4double thePhiStart, 
				const G4double theDeltaPhi, 
				const G4bool thePhiIsOpen, 
				const G4bool isAllBehind )
{
	//
	// Record values
	//
	r[0] = tail->r; z[0] = tail->z;
	r[1] = head->r; z[1] = head->z;
	
	phiIsOpen = thePhiIsOpen;
	if (phiIsOpen) {
		deltaPhi = theDeltaPhi;
		startPhi = thePhiStart;

		//
		// Set phi values to our conventions
		//
		while (deltaPhi < 0.0) deltaPhi += 2.0*M_PI;
		while (startPhi < 0.0) startPhi += 2.0*M_PI;
		
		//
		// Calculate corner coordinates
		//
		corners = new G4ThreeVector[4];
		
		corners[0] = G4ThreeVector( tail->r*cos(startPhi), tail->r*sin(startPhi), tail->z );
		corners[1] = G4ThreeVector( head->r*cos(startPhi), head->r*sin(startPhi), head->z );
		corners[2] = G4ThreeVector( tail->r*cos(startPhi+deltaPhi), tail->r*sin(startPhi+deltaPhi), tail->z );
		corners[3] = G4ThreeVector( head->r*cos(startPhi+deltaPhi), head->r*sin(startPhi+deltaPhi), head->z );
	}
	else {
		deltaPhi = 2*M_PI;
		startPhi = 0.0;
	}
	
	allBehind = isAllBehind;
		
	//
	// Make our intersecting cone
	//
	cone = new G4IntersectingCone( r, z );
	
	//
	// Calculate vectors in r,z space
	//
	rS = r[1]-r[0]; zS = z[1]-z[0];
	length = sqrt( rS*rS + zS*zS);
	rS /= length; zS /= length;
	
	rNorm = +zS;
	zNorm = -rS;
	
	G4double lAdj;
	
	prevRS = r[0]-prevRZ->r;
	prevZS = z[0]-prevRZ->z;
	lAdj = sqrt( prevRS*prevRS + prevZS*prevZS );
	prevRS /= lAdj;
	prevZS /= lAdj;

	rNormEdge[0] = rNorm + prevZS;
	zNormEdge[0] = zNorm - prevRS;
	lAdj = sqrt( rNormEdge[0]*rNormEdge[0] + zNormEdge[0]*zNormEdge[0] );
	rNormEdge[0] /= lAdj;
	zNormEdge[0] /= lAdj;

	nextRS = nextRZ->r-r[1];
	nextZS = nextRZ->z-z[1];
	lAdj = sqrt( nextRS*nextRS + nextZS*nextZS );
	nextRS /= lAdj;
	nextZS /= lAdj;

	rNormEdge[1] = rNorm + nextZS;
	zNormEdge[1] = zNorm - nextRS;
	lAdj = sqrt( rNormEdge[1]*rNormEdge[1] + zNormEdge[1]*zNormEdge[1] );
	rNormEdge[1] /= lAdj;
	zNormEdge[1] /= lAdj;
}


//
// Destructor
//	
G4PolyconeSide::~G4PolyconeSide()
{
	delete cone;
	if (phiIsOpen) delete [] corners;
}


//
// Copy constructor
//
G4PolyconeSide::G4PolyconeSide( const G4PolyconeSide &source )
{
	CopyStuff( source );
}


//
// Assignment operator
//
G4PolyconeSide *G4PolyconeSide::operator=( const G4PolyconeSide &source )
{
	if (this == &source) return this;

	delete cone;
	if (phiIsOpen) delete [] corners;
	
	CopyStuff( source );
	
	return this;
}


//
// CopyStuff
//
void G4PolyconeSide::CopyStuff( const G4PolyconeSide &source )
{
	r[0]		= source.r[0];
	r[1]		= source.r[1];
	z[0]		= source.z[0];
	z[1]		= source.z[1];
	
	startPhi	= source.startPhi;
	deltaPhi	= source.deltaPhi;
	phiIsOpen	= source.phiIsOpen;
	allBehind	= source.allBehind;
	
	cone		= new G4IntersectingCone( *source.cone );
	
	rNorm		= source.rNorm;
	zNorm		= source.zNorm;
	rS		= source.rS;
	zS		= source.zS;
	length		= source.length;
	prevRS		= source.prevRS;
	prevZS		= source.prevZS;
	nextRS		= source.nextRS;
	nextZS		= source.nextZS;
	
	rNormEdge[0] 	= source.rNormEdge[0];
	rNormEdge[1]	= source.rNormEdge[1];
	zNormEdge[0]	= source.zNormEdge[0];
	zNormEdge[1]	= source.zNormEdge[1];
	
	if (phiIsOpen) {
		corners = new G4ThreeVector[4];
		
		corners[0] = source.corners[0];
		corners[1] = source.corners[1];
		corners[2] = source.corners[2];
		corners[3] = source.corners[3];
	}
}


//
// Intersect
//
G4bool G4PolyconeSide::Intersect( const G4ThreeVector &p, const G4ThreeVector &v,	
				  const G4bool outgoing, const G4double surfTolerance,
				  G4double &distance, G4double &distFromSurface,
				  G4ThreeVector &normal, G4bool &isAllBehind )
{
	G4double s1, s2;
	G4double normSign = outgoing ? +1 : -1;
	
	isAllBehind = allBehind;

	//
	// Check for two possible intersections
	//
	G4int nside = cone->LineHitsCone( p, v, &s1, &s2 );
	if (nside == 0) return false;
		
	//
	// Check the first side first, since it is (supposed to be) closest
	//
	G4ThreeVector hit = p + s1*v;
	
	if (PointOnCone( hit, normSign, p, v, normal )) {
		//
		// Good intersection! What about the normal? 
		//
		if (normSign*v.dot(normal) > 0) {
			//
			// We have a valid intersection, but it could very easily
			// be behind the point. To decide if we tolerate this,
			// we have to see if the point p is on the surface near
			// the intersecting point.
			//
			// What does it mean exactly for the point p to be "near"
			// the intersection? It means that if we draw a line from
			// p to the hit, the line remains entirely within the
			// tolerance bounds of the cone. To test this, we can
			// ask if the normal is correct near p.
			//
			G4double pr = p.perp();
			if (pr < DBL_MIN) pr = DBL_MIN;
			G4ThreeVector pNormal( rNorm*p.x()/pr, rNorm*p.y()/pr, zNorm );
			if (normSign*v.dot(pNormal) > 0) {
				//
				// p and intersection in same hemisphere
				//
				G4double distOutside2;
				distFromSurface = -normSign*DistanceAway( p, false, distOutside2 );
				if (distOutside2 < surfTolerance*surfTolerance) {
					if (distFromSurface > -surfTolerance) {
						//
						// We are just inside or away from the
						// surface. Accept *any* value of distance.
						//
						distance = s1;
						return true;
					}
				}
			}
			else 
				distFromSurface = s1;
			
			//
			// Accept positive distances
			//
			if (s1 > 0) {
				distance = s1;
				return true;
			}
		}
	}	
	
	if (nside==1) return false;
	
	//
	// Well, try the second hit
	//	
	hit = p + s2*v;
	
	if (PointOnCone( hit, normSign, p, v, normal )) {
		//
		// Good intersection! What about the normal? 
		//
		if (normSign*v.dot(normal) > 0) {
			G4double pr = p.perp();
			if (pr < DBL_MIN) pr = DBL_MIN;
			G4ThreeVector pNormal( rNorm*p.x()/pr, rNorm*p.y()/pr, zNorm );
			if (normSign*v.dot(pNormal) > 0) {
				G4double distOutside2;
				distFromSurface = -normSign*DistanceAway( p, false, distOutside2 );
				if (distOutside2 < surfTolerance*surfTolerance) {
					if (distFromSurface > -surfTolerance) {
						distance = s2;
						return true;
					}
				}
			}
			else 
				distFromSurface = s2;
			
			if (s2 > 0) {
				distance = s2;
				return true;
			}
		}
	}	

	//
	// Better luck next time
	//
	return false;
}


G4double G4PolyconeSide::Distance( const G4ThreeVector &p, const G4bool outgoing )
{
	G4double normSign = outgoing ? -1 : +1;
	G4double distFrom, distOut2;
	
	//
	// We have two tries for each hemisphere. Try the closest first.
	//
	distFrom = normSign*DistanceAway( p, false, distOut2 );
	if (distFrom > -0.5*kCarTolerance ) {
		//
		// Good answer
		//
		if (distOut2 > 0) 
			return sqrt( distFrom*distFrom + distOut2 );
		else 
			return fabs(distFrom);
	}
	
	//
	// Try second side. 
	//
	distFrom = normSign*DistanceAway( p,  true, distOut2 );
	if (distFrom > -0.5*kCarTolerance) {

		if (distOut2 > 0) 
			return sqrt( distFrom*distFrom + distOut2 );
		else
			return fabs(distFrom);
	}
	
	return kInfinity;
}


//
// Inside
//
EInside G4PolyconeSide::Inside( const G4ThreeVector &p, const G4double tolerance, 
				G4double *bestDistance )
{
	//
	// Check both sides
	//
	G4double distFrom[2], distOut2[2], dist2[2];
	G4double edgeRZnorm[2];
		 
	distFrom[0] =  DistanceAway( p, false, distOut2[0], edgeRZnorm   );
	distFrom[1] =  DistanceAway( p,  true, distOut2[1], edgeRZnorm+1 );
	
	dist2[0] = distFrom[0]*distFrom[0] + distOut2[0];
	dist2[1] = distFrom[1]*distFrom[1] + distOut2[1];
	
	//
	// Who's closest?
	//
	G4int i = fabs(dist2[0]) < fabs(dist2[1]) ? 0 : 1;
	
	*bestDistance = sqrt( dist2[i] );
	
	//
	// Okay then, inside or out?
	//
	if ( (fabs(edgeRZnorm[i]) < tolerance) && (distOut2[i] < tolerance*tolerance) )
		return kSurface;
	else if (edgeRZnorm[i] < 0)
		return kInside;
	else
		return kOutside;
}


//
// Normal
//
G4ThreeVector G4PolyconeSide::Normal( const G4ThreeVector &p,  G4double *bestDistance )
{
	G4ThreeVector dFrom;
	G4double dOut2;
	
	dFrom = DistanceAway( p, false, dOut2 );
	
	*bestDistance = sqrt( dFrom*dFrom + dOut2 );
	
	G4double rad = p.perp();
	return G4ThreeVector( rNorm*p.x()/rad, rNorm*p.y()/rad, zNorm );
}


//
// Extent
//
G4double G4PolyconeSide::Extent( const G4ThreeVector axis )
{
	if (axis.perp2() < DBL_MIN) {
		//
		// Special case
		//
		return axis.z() < 0 ? -cone->ZLo() : cone->ZHi();
	}

	//
	// Is the axis pointing inside our phi gap?
	//
	if (phiIsOpen) {
		G4double phi = axis.phi();
		while( phi < startPhi ) phi += 2*M_PI;
		
		if (phi > deltaPhi+startPhi) {
			//
			// Yeah, looks so. Make four three vectors defining the phi
			// opening
			//
			G4double cosP = cos(startPhi), sinP = sin(startPhi);
			G4ThreeVector a( r[0]*cosP, r[0]*sinP, z[0] );
			G4ThreeVector b( r[1]*cosP, r[1]*sinP, z[1] );
			cosP = cos(startPhi+deltaPhi); sinP = sin(startPhi+deltaPhi);
			G4ThreeVector c( r[0]*cosP, r[0]*sinP, z[0] );
			G4ThreeVector d( r[1]*cosP, r[1]*sinP, z[1] );
			
			G4double ad = axis.dot(a),
				 bd = axis.dot(b),
				 cd = axis.dot(c),
				 dd = axis.dot(d);
			
			if (bd > ad) ad = bd;
			if (cd > ad) ad = cd;
			if (dd > ad) ad = dd;
			
			return ad;
		}
	}

	//
	// Check either end
	//
	G4double aPerp = axis.perp();
	
	G4double a = aPerp*r[0] + axis.z()*z[0];
	G4double b = aPerp*r[1] + axis.z()*z[1];
	
	if (b > a) a = b;
	
	return a;
}



//
// CalculateExtent
//
// See notes in G4VCSGface
//
void G4PolyconeSide::CalculateExtent( const EAxis axis, 
				      const G4VoxelLimits &voxelLimit,
				      const G4AffineTransform &transform,
				      G4SolidExtentList &extentList        )
{
	G4ClippablePolygon polygon;
	
	//
	// Here we will approximate (ala G4Cons) and divide our conical section
	// into segments, like G4Polyhedra. When doing so, the radius
	// is extented far enough such that the segments always lie
	// just outside the surface of the conical section we are
	// approximating.
	//
	
	//
	// Choose phi size of our segment(s) based on constants as
	// defined in meshdefs.hh
	//
	G4int numPhi = (G4int)(deltaPhi/kMeshAngleDefault) + 1;
	if (numPhi < kMinMeshSections) 
		numPhi = kMinMeshSections;
	else if (numPhi > kMaxMeshSections)
		numPhi = kMaxMeshSections;
		
	G4double sigPhi = deltaPhi/numPhi;
	
	//
	// Determine radius factor to keep segments outside
	//
	G4double rFudge = 1.0/cos(0.5*sigPhi);
	
	//
	// Decide which radius to use on each end of the side,
	// and whether a transition mesh is required
	//
	// {r0,z0}  - Beginning of this side
	// {r1,z1}  - Ending of this side
	// {r2,z0}  - Beginning of transition piece connecting previous
	//            side (and ends at beginning of this side)
	//
	// So, order is 2 --> 0 --> 1.
	//                    -------
	//
	// r2 < 0 indicates that no transition piece is required
	//
	G4double r0, r1, r2, z0, z1;
	
	r2 = -1;	// By default: no transition piece
	
	if (rNorm < -DBL_MIN) {
		//
		// This side faces *inward*, and so our mesh has
		// the same radius
		//
		r1 = r[1];
		z1 = z[1];
		z0 = z[0];
		r0 = r[0];
		
		r2 = -1;
		
		if (prevZS > DBL_MIN) {
			//
			// The previous side is facing outwards
			//
			if ( prevRS*zS - prevZS*rS > 0 ) {
				//
				// Transition was convex: build transition piece
				//
				if (r[0] > DBL_MIN) r2 = r[0]*rFudge;
			}
			else {
				//
				// Transition was concave: short this side
				//
				FindLineIntersect( z0, r0, zS, rS,
						   z0, r0*rFudge, prevZS, prevRS*rFudge, z0, r0 );
			}
		}
		
		if ( nextZS > DBL_MIN && (rS*nextZS - zS*nextRS < 0) ) {
			//
			// The next side is facing outwards, forming a 
			// concave transition: short this side
			//
			FindLineIntersect( z1, r1, zS, rS,
					   z1, r1*rFudge, nextZS, nextRS*rFudge, z1, r1 );
		}
	}
	else if (rNorm > DBL_MIN) {
		//
		// This side faces *outward* and is given a boost to
		// it radius
		//
		r0 = r[0]*rFudge;
		z0 = z[0];
		r1 = r[1]*rFudge;
		z1 = z[1];
		
		if (prevZS < -DBL_MIN) {
			//
			// The previous side is facing inwards
			//
			if ( prevRS*zS - prevZS*rS > 0 ) {
				//
				// Transition was convex: build transition piece
				//
				if (r[0] > DBL_MIN) r2 = r[0];
			}
			else {
				//
				// Transition was concave: short this side
				//
				FindLineIntersect( z0, r0, zS, rS*rFudge,
						   z0, r[0], prevZS, prevRS, z0, r0 );
			}
		}
		
		if ( nextZS < -DBL_MIN && (rS*nextZS - zS*nextRS < 0) ) {
			//
			// The next side is facing inwards, forming a 
			// concave transition: short this side
			//
			FindLineIntersect( z1, r1, zS, rS*rFudge,
					   z1, r[1], nextZS, nextRS, z1, r1 );
		}
	}
	else {
		//
		// This side is perpendicular to the z axis (is a disk)
		//
		// Whether or not r0 needs a rFudge factor depends
		// on the normal of the previous edge. Similar with r1
		// and the next edge. No transition piece is required.
		//
		r0 = r[0];
		r1 = r[1];
		z0 = z[0];
		z1 = z[1];
		
		if (prevZS > DBL_MIN) r0 *= rFudge;
		if (nextZS > DBL_MIN) r1 *= rFudge;
	}
	
	//
	// Loop
	//
	G4double phi = startPhi, 
	         cosPhi = cos(phi), 
		 sinPhi = sin(phi);
	
	G4ThreeVector v0( r0*cosPhi, r0*sinPhi, z0 ),
		      v1( r1*cosPhi, r1*sinPhi, z1 ),
		      v2, w0, w1, w2;
	transform.ApplyPointTransform( v0 );
	transform.ApplyPointTransform( v1 );
	
	if (r2 >= 0) {
		v2 = G4ThreeVector( r2*cosPhi, r2*sinPhi, z0 );
		transform.ApplyPointTransform( v2 );
	}

	do {
		G4double min, max;
	
		phi += sigPhi;
		if (numPhi == 1) phi = startPhi+deltaPhi;	// Try to avoid roundoff
	        cosPhi = cos(phi), 
		sinPhi = sin(phi);
		
		w0 = G4ThreeVector( r0*cosPhi, r0*sinPhi, z0 );
		w1 = G4ThreeVector( r1*cosPhi, r1*sinPhi, z1 );
		transform.ApplyPointTransform( w0 );
		transform.ApplyPointTransform( w1 );
		
		G4ThreeVector deltaV = r0 > r1 ? w0-v0 : w1-v1;
		
		//
		// Build polygon, taking special care to keep the vertices
		// in order
		//
		polygon.ClearAllVertices();

		polygon.AddVertexInOrder( v0 );
		polygon.AddVertexInOrder( v1 );
		polygon.AddVertexInOrder( w1 );
		polygon.AddVertexInOrder( w0 );

		//
		// Get extent
		//
		if (polygon.PartialClip( voxelLimit, axis )) {
			//
			// Get dot product of normal with target axis
			//
			polygon.SetNormal( deltaV.cross(v1-v0).unit() );
			
			extentList.AddSurface( polygon );
		}
		
		if (r2 >= 0) {
			//
			// Repeat, for transition piece
			//
			w2 = G4ThreeVector( r2*cosPhi, r2*sinPhi, z0 );
			transform.ApplyPointTransform( w2 );

			polygon.ClearAllVertices();

			polygon.AddVertexInOrder( v2 );
			polygon.AddVertexInOrder( v0 );
			polygon.AddVertexInOrder( w0 );
			polygon.AddVertexInOrder( w2 );

			if (polygon.PartialClip( voxelLimit, axis )) {
				polygon.SetNormal( deltaV.cross(v0-v2).unit() );
				
				extentList.AddSurface( polygon );
			}
			
			v2 = w2;
		}
		
		//
		// Next vertex
		//		
		v0 = w0;
		v1 = w1;
	} while( --numPhi > 0 );
	
	//
	// We are almost done. But, it is important that we leave no
	// gaps in the surface of our solid. By using rFudge, however,
	// we've done exactly that, if we have a phi segment. 
	// Add two additional faces if necessary
	//
	if (phiIsOpen && rNorm > DBL_MIN) {
		G4double min, max;
		
		G4double cosPhi = cos(startPhi),
			 sinPhi = sin(startPhi);

		G4ThreeVector a0( r[0]*cosPhi, r[0]*sinPhi, z[0] ),
			      a1( r[1]*cosPhi, r[1]*sinPhi, z[1] ),
			      b0( r0*cosPhi, r0*sinPhi, z[0] ),
			      b1( r1*cosPhi, r1*sinPhi, z[1] );
	
		transform.ApplyPointTransform( a0 );
		transform.ApplyPointTransform( a1 );
		transform.ApplyPointTransform( b0 );
		transform.ApplyPointTransform( b1 );

		polygon.ClearAllVertices();

		polygon.AddVertexInOrder( a0 );
		polygon.AddVertexInOrder( a1 );
		polygon.AddVertexInOrder( b0 );
		polygon.AddVertexInOrder( b1 );
		
		if (polygon.PartialClip( voxelLimit , axis)) {
			G4ThreeVector normal( sinPhi, -cosPhi, 0 );
			polygon.SetNormal( transform.TransformAxis( normal ) );
				
			extentList.AddSurface( polygon );
		}
		
		cosPhi = cos(startPhi+deltaPhi);
		sinPhi = sin(startPhi+deltaPhi);
		
		a0 = G4ThreeVector( r[0]*cosPhi, r[0]*sinPhi, z[0] ),
		a1 = G4ThreeVector( r[1]*cosPhi, r[1]*sinPhi, z[1] ),
		b0 = G4ThreeVector( r0*cosPhi, r0*sinPhi, z[0] ),
		b1 = G4ThreeVector( r1*cosPhi, r1*sinPhi, z[1] );
		transform.ApplyPointTransform( a0 );
		transform.ApplyPointTransform( a1 );
		transform.ApplyPointTransform( b0 );
		transform.ApplyPointTransform( b1 );

		polygon.ClearAllVertices();

		polygon.AddVertexInOrder( a0 );
		polygon.AddVertexInOrder( a1 );
		polygon.AddVertexInOrder( b0 );
		polygon.AddVertexInOrder( b1 );
		
		if (polygon.PartialClip( voxelLimit, axis )) {
			G4ThreeVector normal( -sinPhi, cosPhi, 0 );
			polygon.SetNormal( transform.TransformAxis( normal ) );
				
			extentList.AddSurface( polygon );
		}
	}
		
	return;
}


//
// -------------------------------------------------------

//
// DistanceAway
//
// Calculate distance of a point from our conical surface, including the effect
// of any phi segmentation
//
// Arguments:
//	p		- (in) Point to check
//	opposite	- (in) If true, check opposite hemisphere (see below)
//	distOutside	- (out) Additional distance outside the edges of the
//			 	surface
//	edgeRZnorm	- (out) if negative, point is inside
//	return value = distance from the conical plane, if extrapolated beyond edges,
//		       signed by whether the point is in inside or outside the shape
//
// Notes:
//	* There are two answers, depending on which hemisphere is considered.
//
G4double G4PolyconeSide::DistanceAway( const G4ThreeVector &p, const G4bool opposite,
				       G4double &distOutside2, G4double *edgeRZnorm  )
{
	//
	// Convert our point to r and z
	//
	G4double rx = p.perp(), zx = p.z();
	
	//
	// Change sign of r if opposite says we should
	//
	if (opposite) rx = -rx;
	
	//
	// Calculate return value
	//
	G4double deltaR  = rx - r[0], deltaZ = zx - z[0];
	G4double answer = deltaR*rNorm + deltaZ*zNorm;
	
	//
	// Are we off the surface in r,z space?
	//
	G4double s = deltaR*rS + deltaZ*zS;
	if (s < 0) {
		distOutside2 = s*s;
		if (edgeRZnorm) *edgeRZnorm = deltaR*rNormEdge[0] + deltaZ*zNormEdge[0];
	}
	else if (s > length) {
		distOutside2 = sqr( s-length );
		if (edgeRZnorm) {
			G4double deltaR  = rx - r[1], deltaZ = zx - z[1];
			*edgeRZnorm = deltaR*rNormEdge[1] + deltaZ*zNormEdge[1];
		}
	}
	else {
		distOutside2 = 0;
		if (edgeRZnorm) *edgeRZnorm = answer;
	}

	if (phiIsOpen) {
		//
		// Finally, check phi
		//
		G4double phi = p.phi();
		while( phi < startPhi ) phi += 2*M_PI;
		
		if (phi > startPhi+deltaPhi) {
			//
			// Oops. Are we closer to the start phi or end phi?
			//
			G4double d1 = phi-startPhi-deltaPhi;
			while( phi > startPhi ) phi -= 2*M_PI;
			G4double d2 = startPhi-phi;
			
			if (d2 < d1) d1 = d2;
			
			//
			// Add result to our distance
			//
			G4double dist = d1*rx;
			
			distOutside2 += dist*dist;
			if (edgeRZnorm) *edgeRZnorm = fabs(dist);
		}
	}

	return answer;
}


//
// PointOnCone
//
// Decide if a point is on a cone and return normal if it is
//
G4bool G4PolyconeSide::PointOnCone( const G4ThreeVector &hit, const G4double normSign,
				    const G4ThreeVector &p, const G4ThreeVector &v,
			 	    G4ThreeVector &normal )
{
	G4double rx = hit.perp();
	//
	// Check radial/z extent, as appropriate
	//
	if (!cone->HitOn( rx, hit.z() )) return false;
	
	if (phiIsOpen) {
		G4double phiTolerant = 2.0*kCarTolerance/(rx+kCarTolerance);
		//
		// Check phi segment. Here we have to be careful
		// to use the standard method consistent with
		// PolyPhiFace. See PolyPhiFace::InsideEdgesExact
		//
		G4double phi = hit.phi();
		while( phi < startPhi-phiTolerant ) phi += 2*M_PI;
		
		if (phi > startPhi+deltaPhi+phiTolerant) return false;
		
		if (phi > startPhi+deltaPhi-phiTolerant) {
			//
			// Exact treatment
			//
			G4ThreeVector qx = p + v;
			G4ThreeVector qa = qx - corners[2],
				      qb = qx - corners[3];
			G4ThreeVector qacb = qa.cross(qb);
			
			if (normSign*qacb.dot(v) < 0) return false;
		}
		else if (phi < phiTolerant) {
			G4ThreeVector qx = p + v;
			G4ThreeVector qa = qx - corners[1],
				      qb = qx - corners[0];
			G4ThreeVector qacb = qa.cross(qb);
			
			if (normSign*qacb.dot(v) < 0) return false;
		}
	}
	
	//
	// We have a good hit! Calculate normal
	//
	if (rx < DBL_MIN) 
		normal = G4ThreeVector( 0, 0, zNorm < 0 ? -1 : 1 );
	else
		normal = G4ThreeVector( rNorm*hit.x()/rx, rNorm*hit.y()/rx, zNorm );
	return true;
}


//
// FindLineIntersect
//
// Decide the point at which two 2-dimensional lines intersect
//
// Equation of line: x = x1 + s*tx1
//                   y = y1 + s*ty1
//
// It is assumed that the lines are *not* parallel
//
void G4PolyconeSide::FindLineIntersect( const G4double x1, const G4double y1,
					const G4double tx1, const G4double ty1,
				        const G4double x2, const G4double y2,
					const G4double tx2, const G4double ty2,
					G4double &x, G4double &y )
{
	//
	// The solution is a simple linear equation
	//
	G4double deter = tx1*ty2 - tx2*ty1;
	
	G4double s1 = ((x2-x1)*ty2 - tx2*(y2-y1))/deter;
	G4double s2 = ((x2-x1)*ty1 - tx1*(y2-y1))/deter;

	//
	// We want the answer to not depend on which order the
	// lines were specified. Take average.
	//
	x = 0.5*( x1+s1*tx1 + x2+s2*tx2 );
	y = 0.5*( y1+s1*ty1 + y2+s2*ty2 );
}
