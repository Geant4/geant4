//
// G4PolyconeSide.cc
//
// Implemenation of the face representing one conical side of a polycone
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4PolyconeSide.hh"
#include "G4IntersectingCone.hh"
#include "G4ClippablePolygon.hh"
#include "G4AffineTransform.hh"
#include "meshdefs.hh"

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
				const G4bool thePhiIsOpen )
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
	}
	else {
		deltaPhi = 2*M_PI;
		startPhi = 0.0;
	}
		
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
	
	if (r[0] < 1/kInfinity) {
		//
		// This segment begins at R=0 
		//
		prevRNorm = 0;
		rNormEdge[0] = rNorm;
		zNormEdge[0] = zNorm;
	}
	else {
		G4double rAdj = r[0]-prevRZ->r, zAdj = z[0]-prevRZ->z;
		G4double lAdj = sqrt( rAdj*rAdj + zAdj*zAdj );
		rAdj /= lAdj;
		zAdj /= lAdj;

		prevRNorm = zAdj;

		rNormEdge[0] = rNorm + zAdj;
		zNormEdge[0] = zNorm - rAdj;
		lAdj = sqrt( rNormEdge[0]*rNormEdge[0] + zNormEdge[0]*zNormEdge[0] );
		rNormEdge[0] /= lAdj;
		zNormEdge[0] /= lAdj;
	}
	
	if (r[1] < 1/kInfinity) {
		//
		// This segment ends at R=0
		//
		rNormEdge[1] = rNorm;
		zNormEdge[1] = zNorm;
	}
	else {
		G4double rAdj = nextRZ->r-r[1], zAdj = nextRZ->z-z[1];
		G4double lAdj = sqrt( rAdj*rAdj + zAdj*zAdj );
		rAdj /= lAdj;
		zAdj /= lAdj;

		rNormEdge[1] = rNorm + zAdj;
		zNormEdge[1] = zNorm - rAdj;
		lAdj = sqrt( rNormEdge[1]*rNormEdge[1] + zNormEdge[1]*zNormEdge[1] );
		rNormEdge[1] /= lAdj;
		zNormEdge[1] /= lAdj;
	}
}


//
// Destructor
//	
G4PolyconeSide::~G4PolyconeSide()
{
	delete cone;
}


G4bool G4PolyconeSide::Intersect( const G4ThreeVector &p, const G4ThreeVector &v,	
				  const G4bool outgoing, const G4double surfTolerance,
				  G4double &distance, G4double &distFromSurface,
				  G4ThreeVector &normal, G4bool &allBehind )
{
	G4double s1, s2;
	G4double normSign = outgoing ? +1 : -1;
	
	allBehind = true;

	//
	// Check for two possible intersections
	//
	G4int nside = cone->LineHitsCone( p, v, &s1, &s2 );
	if (nside == 0) return false;
	
	//
	// Check the first side first, since it is (supposed to be) closest
	//
	G4ThreeVector hit = p + s1*v;
	
	if (PointOnCone( hit, normal )) {
		//
		// Good intersection! What about the normal? 
		//
		G4double vdotN = v.dot(normal);
		if (normSign*vdotN > 0) {
			//
			// Good normal! Check distance
			//
			// We want to apply an addition surface tolerance
			// if the point p is on the surface. Is it?
			//
			G4bool opposite = (p.x()*hit.x()+p.y()*hit.y() < 0);
			G4double distOutside2;
			G4double notUsed = -normSign*DistanceAway( p, opposite, distOutside2, 0 );
					
			//
			// The distance from the surface is defined along the
			// normal at the intersection point.
			//
			// I think.
			//	
			distFromSurface = normSign*s1*vdotN;
		
			//
			// Apply tolerance, but only if the point is outside 
			// the edges of the cone
			//
			if (distFromSurface > (distOutside2 > 0 ? 0 : -surfTolerance)) {
				//
				// Good intersection. Return now, since it is the closest.
				//
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
	
	if (PointOnCone( hit, normal )) {
		//
		// Good intersection! What about the normal? 
		//
		G4double vdotN = v.dot(normal);
		if (normSign*vdotN > 0) {
			//
			// Good normal! Check distance
			//
			// We want to apply an addition surface tolerance
			// if the point p is on the surface. Is it?
			//
			G4bool opposite = (p.x()*hit.x()+p.y()*hit.y() < 0);
			G4double distOutside2;
			G4double notUsed = -normSign*DistanceAway( p, opposite, distOutside2, 0 );
					
			//
			// The distance from the surface is defined along the
			// normal at the intersection point.
			//
			// I think.
			//	
			distFromSurface = normSign*s2*vdotN;
		
			//
			// Apply tolerance, but only if the point is outside 
			// the edges of the cone
			//
			if (distFromSurface > (distOutside2 > 0 ? 0 : -surfTolerance)) {
				//
				// Good intersection. Return now, since it is the closest.
				//
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
	distFrom = DistanceAway( p, false, distOut2, 0 );
	if (distFrom*normSign > 0) {
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
	distFrom = DistanceAway( p,  true, distOut2, 0 );
	if (distFrom*normSign > 0) {

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
		 
	distFrom[0] =  DistanceAway( p, false, distOut2[0], edgeRZnorm );
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
	
	dFrom = DistanceAway( p, false, dOut2, 0 );
	
	*bestDistance = sqrt( dFrom*dFrom + dOut2 );
	
	G4double rad = p.perp();
	return G4ThreeVector( rNorm*p.x()/rad, rNorm*p.y()/rad, zNorm );
}


//
// Extent
//
G4double G4PolyconeSide::Extent( const G4ThreeVector axis )
{
	if (axis.perp2() < 1.0/kInfinity) {
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
				      G4double &min, G4double &max        )
{
	G4ClippablePolygon polygon;
	
	//
	// Here we will cheat (ala G4Cons) and divide our conical section
	// into segments, like G4Polyhedra. When doing so, the radius
	// is extented far enough such that the segments always lie
	// just outside the surface of the conical section we are
	// approximating.
	//
	
	//
	// Choose phi size of our segment(s) based on constants as
	// defined in meshdefs.hh
	//
	G4int numPhi = deltaPhi/kMeshAngleDefault + 1;
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
	G4double r0, r1, r2;
	
	if (rNorm < -1/kInfinity) {
		//
		// This side faces *inward*
		//
		r0 = r[0];
		r1 = r[1];
		
		//
		// A transition is required if the previous side
		// faced outward
		//
		r2 = (prevRNorm > 1/kInfinity) ? r0*rFudge : -1;
	}
	else if (rNorm > 1/kInfinity) {
		//
		// This side faces *outward*
		//
		r0 = r[0]*rFudge;
		r1 = r[1]*rFudge;
		
		//
		// A transition is required if the previous side
		// faced inward
		//
		r2 = (prevRNorm < -1/kInfinity) ? r[0] : -1;
	}
	else {
		//
		// This side is perpendicular to the z axis (is a disk)
		//
		r0 = r[0];
		r1 = r[1];
		if (r0 > r1) r0 *= rFudge; else r1 *= rFudge;
	}
	
	//
	// Loop
	//
	G4double phi = startPhi, 
	         cosPhi = cos(phi), 
		 sinPhi = sin(phi);
	
	G4ThreeVector v0( r0*cosPhi, r0*sinPhi, z[0] ),
		      v1( r1*cosPhi, r1*sinPhi, z[1] ),
		      v2( r2*cosPhi, r2*sinPhi, z[0] ),
		      w0, w1, w2;
	transform.ApplyPointTransform( v0 );
	transform.ApplyPointTransform( v1 );
	transform.ApplyPointTransform( v2 );

	do {
		phi += sigPhi;
	        cosPhi = cos(phi), 
		sinPhi = sin(phi);
		
		w0 = G4ThreeVector( r0*cosPhi, r0*sinPhi, z[0] );
		w1 = G4ThreeVector( r1*cosPhi, r1*sinPhi, z[1] );
		transform.ApplyPointTransform( w0 );
		transform.ApplyPointTransform( w1 );
		
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
		polygon.Clip( voxelLimit );
		polygon.GetExtent( axis, min, max );
		
		if (r2 >= 0) {
			//
			// Repeat, for transition piece
			//
			w2 = G4ThreeVector( r2*cosPhi, r2*sinPhi, z[0] );
			transform.ApplyPointTransform( w2 );

			polygon.ClearAllVertices();
		
			polygon.AddVertexInOrder( v2 );
			polygon.AddVertexInOrder( v0 );
			polygon.AddVertexInOrder( w0 );
			polygon.AddVertexInOrder( w2 );
			
			polygon.Clip( voxelLimit );
			polygon.GetExtent( axis, min, max );
			
			v2 = w2;
		}
		
		//
		// Next vertex
		//		
		v0 = w0;
		v1 = w1;
	} while( --numPhi > 0 );
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
//	edgeNorm	- (out) Edge Status (belowRZ, aboveRZ, inRZ)
//	return value = distance from the conical plane, if extrapolated beyond edges,
//		       signed by whether the point is in inside or outside the shape
//
// Notes:
//	* There are two answers, depending on which hemisphere is considered.
//
G4double G4PolyconeSide::DistanceAway( const G4ThreeVector &p, const G4bool opposite,
				       G4double &distOutside2, G4double *edgeRZnorm )
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
		if (edgeRZnorm) *edgeRZnorm = deltaR*rNormEdge[1] + deltaZ*zNormEdge[1];
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
			distOutside2 += d1*d1*p.perp2();
		}
	}

	return answer;
}


//
// PointOnCone
//
// Decide if a point is on a cone and return normal if it is
//
G4bool G4PolyconeSide::PointOnCone( const G4ThreeVector &p, G4ThreeVector &normal )
{
	G4double rx = p.perp();
	//
	// Check radial/z extent, as appropriate
	//
	if (!cone->HitOn( rx, p.z() )) return false;
	
	if (phiIsOpen) {
		//
		// Check phi segment
		//
		G4double phi = p.phi();
		while( phi < startPhi ) phi += 2*M_PI;
		
		if (phi > startPhi+deltaPhi) return false;
	}
	
	//
	// We have a good hit! Calculate normal
	//
	if (rx<0) rx = p.perp();
	
	if (rx < -1.0/kInfinity) 
		normal = G4ThreeVector( 0, 0, zNorm < 0 ? -1 : 1 );
	else
		normal = G4ThreeVector( rNorm*p.x()/rx, rNorm*p.y()/rx, zNorm );
	return true;
}
