//
// G4PolyPhiFace.cc
//
// Implementation of the face that bounds a polycone or polyhedra at
// its phi opening.
//
// ----------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "G4PolyPhiFace.hh"
#include "G4ClippablePolygon.hh"
#include "G4ReduciblePolygon.hh"
#include "G4AffineTransform.hh"

//
// Constructor
//
// Points r,z should be supplied in clockwise order in r,z. For example:
//
//                [1]---------[2]         ^ R
//                 |           |          |
//                 |           |          +--> z
//                [0]---------[3]
//
G4PolyPhiFace::G4PolyPhiFace( const G4ReduciblePolygon *rz, const G4double phi, 
			      const G4double deltaPhi, const G4double phiOther )
{
	numEdges = rz->NumVertices();
	
	rMin = rz->Amin();
	rMax = rz->Amax();
	zMin = rz->Bmin();
	zMax = rz->Bmax();

	//
	// Is this the "starting" phi edge of the two?
	//
	G4bool start = (phiOther > phi);
	
	//
	// Build radial vector
	//
	radial = G4ThreeVector( cos(phi), sin(phi), 0.0 );
	
	//
	// Build normal
	//
	G4double zSign = start ? 1 : -1;
	normal = G4ThreeVector( zSign*radial.y(), -zSign*radial.x(), 0 );

        //
        // Allocate corners
        //
        corners = new G4PolyPhiFaceVertex[numEdges];

        //
        // Fill them
        //
	G4ReduciblePolygonIterator iterRZ(rz);
	
	G4PolyPhiFaceVertex *corn = corners;
	iterRZ.Begin();
	do {
		corn->r = iterRZ.GetA();
		corn->z = iterRZ.GetB();
	} while( ++corn, iterRZ.Next() );

        //
        // Allocate edges
        //
        edges = new G4PolyPhiFaceEdge[numEdges];

        //
        // Fill them
        //
	G4double midPhi = phi + (start ? +0.5 : -0.5)*deltaPhi;
	G4double cosMid = cos(midPhi), 
		 sinMid = sin(midPhi);
	G4double rFact = cos(0.5*deltaPhi);
	G4ThreeVector sideNorm;

        G4PolyPhiFaceVertex *prev = corners+numEdges-1,
                            *here = corners;
        G4PolyPhiFaceEdge   *edge = edges;
        do {
                edge->v0 = prev;
                edge->v1 = here;

                G4double dr = here->r - prev->r,
                         dz = here->z - prev->z;
                         
                edge->length = sqrt( dr*dr + dz*dz );

                edge->tr = dr/edge->length;
                edge->tz = dz/edge->length;
		
		if ((here->r < 1/kInfinity) && (prev->r < 1/kInfinity)) {
			//
			// Sigh! Always exceptions!
			// This edge runs at r==0, so its adjoing surface is not a
			// PolyconeSide or PolyhedraSide, but the opposite PolyPhiFace.
			//
			G4double zSignOther = start ? -1 : 1;
			sideNorm = G4ThreeVector(  zSignOther*sin(phiOther), 
					          -zSignOther*cos(phiOther), 0 );
		}
		else {
			sideNorm = G4ThreeVector( dz*cosMid, dz*sinMid, -dr*rFact );
			sideNorm = sideNorm.unit();
		}
		sideNorm += normal;
		
		//
		// --- Just for now --- Leave these normals unnormalized
		//
		edge->norm3D = sideNorm;
        } while( edge++, prev=here, ++here < corners+numEdges );

        //
        // Go back an fill in corner "normals", which are just the
        // average of the normals of the ajoining edges
        //
        G4PolyPhiFaceEdge *prevEdge = edges+numEdges-1;
        edge = edges;
        do {
		G4double rPart = prevEdge->tr + edge->tr;
		G4double zPart = prevEdge->tz + edge->tz;
		G4double norm = sqrt( rPart*rPart + zPart*zPart );
                edge->v0->rNorm = +zPart/norm;
                edge->v0->zNorm = -rPart/norm;

		G4ThreeVector norm3D = prevEdge->norm3D + edge->norm3D - normal;

		if (edge->v0->r < 1/kInfinity &&
		    edge->v1->r > 1/kInfinity &&
		    prevEdge->v0->r > 1/kInfinity ) {
		 	//
			// Oh boy, it doesn't get any more obscure than this!
			// What we have is a corner point at r==0, without
			// any adjacent edges along r==0. Example:
			//
			//                   \****|
			//                    \***|         R
			//                     \**|         ^
			//                      \*|         |
			//                       \|         |
			//  ......(R=0)...........+........ o----> Z
			//
			// Implication? This is a 3D corner that has four
			// adjoining faces, the only possible (phi segmented) 
			// Polycone/hedra corner
			// that has more than three. It is quite a special
			// case, which I will deal with here.
			//
			// (BTW: in such interesting shapes like these, the thickness
			// of the shape on the z axis is zero. The interesting question
			// is: will a particle traveling exactly along the z axis
			// intersect the shape? The reasonable answer: yes. This is
			// dealt with in G4PolyconeSide and G4PolyhedraSide.)
			//
			G4double zSignOther = start ? -1 : 1;
			G4ThreeVector normalOther = G4ThreeVector(  zSignOther*sin(phiOther), 
								   -zSignOther*cos(phiOther), 0 );
								   
			norm3D = prevEdge->norm3D + edge->norm3D + 2*normalOther;

			G4double midPhi = phiOther + (start ? -0.5 : +0.5)*deltaPhi;
			G4double cosMid = cos(midPhi), 
				 sinMid = sin(midPhi);
			G4ThreeVector sideNorm1 = G4ThreeVector( edge->tz*cosMid, 
						 		 edge->tz*sinMid, -edge->tr*rFact );
			G4ThreeVector sideNorm2 = G4ThreeVector( prevEdge->tz*cosMid, 
						 		 prevEdge->tz*sinMid, -prevEdge->tr*rFact );

			norm3D += sideNorm1.unit() + sideNorm2.unit();
		}
		else {
			//
			// Corner normal should be average of normals of connecting edges,
			// or, equivalently, the average of all connecting faces.
			//
			// prevEdge->norm3D = normal + side1.normal = A
			// edge->norm3D     = normal + side2.normal = B
			// A + B - normal = normal + side1.normal + side2.normal
			//
			// This trick only works, though, if the edge->norm3D and 
			// prevEdge->norm3D are left unnormalized!
			//
		
			norm3D = prevEdge->norm3D + edge->norm3D - normal;
		}

		//
		// Normalize the resulting vector
		//
		edge->v0->norm3D = norm3D.unit();
        } while(  prevEdge=edge, ++edge < edges+numEdges );
	
	//
	// Okay, go back and normalize the edge normal vectors
	//
	edge = edges;
	do {
		edge->norm3D = edge->norm3D.unit();
	} while( ++edge < edges+numEdges );

	//
	// Build point on surface
	//
	G4double rAve = 0.5*(rMax-rMin),
		 zAve = 0.5*(zMax-zMin);
	surface = G4ThreeVector( rAve*radial.x(), rAve*radial.y(), zAve );
}


//
// Destructor
//
G4PolyPhiFace::~G4PolyPhiFace()
{
	delete [] edges;
}


//
// Intersect
//
G4bool G4PolyPhiFace::Intersect( const G4ThreeVector &p, const G4ThreeVector &v,
				 const G4bool outgoing, const G4double surfTolerance,
				 G4double &distance, G4double &distFromSurface,
				 G4ThreeVector &aNormal, G4bool &allBehind )
{
	G4double normSign = outgoing ? +1 : -1;
	
	//
	// These don't change
	//
	allBehind = true;
	aNormal = normal;

	//
	// Correct normal? Here we have straight sides, and can safely ignore
	// intersections where the dot product with the normal is zero.
	//
	G4double dotProd = normSign*normal.dot(v);
	
	if (dotProd <= 0) return false;

	//
	// Calculate distance to surface. If the side is too far
	// behind the point, we must reject it.
	//
	G4ThreeVector ps = p - surface;
	distFromSurface = -normSign*ps.dot(normal);
		
	if (distFromSurface < -surfTolerance) return false;

	//
	// Calculate precise distance to intersection with the side
	// (along the trajectory, not normal to the surface)
	//
	distance = distFromSurface/dotProd;

	//
	// Calculate intersection point in r,z
	//
	G4ThreeVector ip = p + distance*v;
	
	G4double r = radial.dot(ip);
	
	//
	// And is it inside the r/z extent?
	//
	return InsideEdges( r, ip.z() );
}


//
// Distance
//
G4double G4PolyPhiFace::Distance( const G4ThreeVector &p, const G4bool outgoing )
{
	G4double normSign = outgoing ? +1 : -1;
	//
	// Correct normal? 
	//
	G4ThreeVector ps = p - surface;
	G4double distPhi = -normSign*normal.dot(ps);
	
	if (distPhi <= 0) return kInfinity;
	
	//
	// Calculate projected point in r,z
	//
	G4double r = radial.dot(p);
	
	//
	// Are we inside the face?
	//
	G4double distRZ2;
	
	if (InsideEdges( r, p.z(), &distRZ2, 0 )) {
		//
		// Yup, answer is just distPhi
		//
		return distPhi;
	}
	else {
		//
		// Nope. Penalize by distance out
		//
		return sqrt( distPhi*distPhi + distRZ2 );
	}
}	
	

//
// Inside
//
EInside G4PolyPhiFace::Inside( const G4ThreeVector &p, const G4double tolerance, 
			       G4double *bestDistance )
{
	//
	// Get distance along phi, which if negative means the point
	// is nominally inside the shape.
	//
	G4ThreeVector ps = p - surface;
	G4double distPhi = normal.dot(ps);
	
	//
	// Calculate projected point in r,z
	//
	G4double r = radial.dot(p);
	
	//
	// Are we inside the face?
	//
	G4double distRZ2;
	G4PolyPhiFaceVertex *base3Dnorm;
	G4ThreeVector	    *head3Dnorm;
	G4bool wereIn = InsideEdges( r, p.z(), &distRZ2, &base3Dnorm, &head3Dnorm );
	
	if (wereIn) {
		//
		// Looks like we're inside. Distance is distance in phi.
		//
		*bestDistance = fabs(distPhi);
	}
	else {
		//
		// We're outside the extent of the face,
		// so the distance is penalized by distance from edges in RZ
		//
		*bestDistance = sqrt( distPhi*distPhi + distRZ2 );
	}
	
	//
	// Can we be on the surface? Yes, but only if we're inside, or
	// close to inside by tolerance
	//	
	if (wereIn || distRZ2 < tolerance*tolerance ) {
		//
		// Yup, answer depends on distPhi, and we can use tolerance
		// to decide if we are on the surface
		//
		if (distPhi < -tolerance) return kInside;
		if (distPhi <  tolerance) return kSurface;
		return kOutside;
	}
	else {
		//
		// Nope. we can only be in or out, and we must
		// used the edge normal to decide
		//
		G4ThreeVector cc( base3Dnorm->r*radial.x(),
				  base3Dnorm->r*radial.y(),
				  base3Dnorm->z );
		cc = p - cc;
		return head3Dnorm->dot(cc) < 0 ? kInside : kOutside;
	}
}	


//
// Normal
//
// This virtual member is simple for our planer shape, which has only one normal
//
G4ThreeVector G4PolyPhiFace::Normal( const G4ThreeVector &p,  G4double *bestDistance )
{
	//
	// Get distance along phi, which if negative means the point
	// is nominally inside the shape.
	//
	G4double distPhi = normal.dot(p);

	//
	// Calculate projected point in r,z
	//
	G4double r = radial.dot(p);
	
	//
	// Are we inside the face?
	//
	G4double distRZ2;
	
	if (InsideEdges( r, p.z(), &distRZ2, 0 )) {
		//
		// Yup, answer is just distPhi
		//
		*bestDistance = fabs(distPhi);
	}
	else {
		//
		// Nope. Penalize by distance out
		//
		*bestDistance = sqrt( distPhi*distPhi + distRZ2 );
	}
	
	return normal;
}


//
// Extent
//
// This actually isn't needed by polycone or polyhedra...
//
G4double G4PolyPhiFace::Extent( const G4ThreeVector axis )
{
	G4double max = -kInfinity;
	
	G4PolyPhiFaceVertex *corner = corners;
	do {
		G4double here = axis.x()*corner->r*radial.x()
			      + axis.y()*corner->r*radial.y()
			      + axis.z()*corner->z;
		if (here > max) max = here;
	} while( ++corner < corners + numEdges );
	
	return max;
}	


//
// CalculateExtent
//
// See notes in G4VCSGface
//
void G4PolyPhiFace::CalculateExtent( const EAxis axis, 
				     const G4VoxelLimits &voxelLimit,
				     const G4AffineTransform &transform,
				     G4double &min, G4double &max        )
{
	//
	// Construct a (sometimes big) clippable polygon, 
	//
	// Perform the necessary transformations while doing so
	//
	G4ClippablePolygon polygon;
	
	G4PolyPhiFaceVertex *corner = corners;
	do {
		G4ThreeVector point( 0, 0, corner->z );
		point += radial*corner->r;
		
		polygon.AddVertexInOrder( transform.TransformPoint( point ) );
	} while( ++corner < corners + numEdges );
	
	//
	// Clip away
	//
	polygon.Clip( voxelLimit );
	
	//
	// Get extent
	//
	polygon.GetExtent( axis, min, max );
}


//
//-------------------------------------------------------
	
//
// InsideEdges (don't care aboud distance)
//
// Decide if the point in r,z is inside the edges of our face
//
// This routine can be made a zillion times quicker by implementing
// better code, for example:
//
//    int pnpoly(int npol, float *xp, float *yp, float x, float y)
//    {
//      int i, j, c = 0;
//      for (i = 0, j = npol-1; i < npol; j = i++) {
//        if ((((yp[i]<=y) && (y<yp[j])) ||
//             ((yp[j]<=y) && (y<yp[i]))) &&
//            (x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
//
//          c = !c;
//      }
//      return c;
//    }
//
// See "Point in Polyon Strategies", Eric Haines [Graphic Gems IV]  pp. 24-46
//
// My algorithm below is rather unique, but is based on code needed to
// calculate the distance to the shape. I left it in here because ...
// well ... to test it better.
//
G4bool G4PolyPhiFace::InsideEdges( const G4double r, const G4double z )
{
	//
	// Quick check of extent
	//
	if ( r < rMin || r > rMax ) return false;
	if ( z < zMin || z > zMax ) return false;
	
	//
	// More thorough check
	//
	G4double notUsed;
	
	return InsideEdges( r, z, &notUsed, 0 );
}


//
// InsideEdges (care about distance)
//
// Decide if the point in r,z is inside the edges of our face
//
G4bool G4PolyPhiFace::InsideEdges( const G4double r, const G4double z,
				   G4double *bestDist2, 
				   G4PolyPhiFaceVertex **base3Dnorm, 
				   G4ThreeVector **head3Dnorm )
{
	G4double bestDistance2 = kInfinity;
	G4bool	 answer;
	
	G4PolyPhiFaceEdge *edge = edges;
	do {
                G4PolyPhiFaceVertex *testMe;
                //
                // Get distance perpendicular to the edge
                //
                G4double dr = (r-edge->v0->r), dz = (z-edge->v0->z);

                G4double distOut = dr*edge->tz - dz*edge->tr;
                G4double distance2 = distOut*distOut;
                if (distance2 > bestDistance2) continue;        // No hope!

                //
                // Check to see if normal intersects edge within the edge's boundary
                //
                G4double s = dr*edge->tr + dz*edge->tz;

		//
		// If it doesn't, penalize distance2 appropriately
		//
		if (s < 0) {
			distance2 += s*s;
                        testMe = edge->v0;
		}
		else if (s > edge->length) {
			G4double s2 = s-edge->length;
			distance2 += s2*s2;
                        testMe = edge->v1;
		}
                else {
                        testMe = 0;
                }
		
		//
		// Closest edge so far?
		//
		if (distance2 < bestDistance2) {
			bestDistance2 = distance2;
                        if (testMe) {
                                G4double distNorm = dr*testMe->rNorm + dz*testMe->zNorm;
                                answer = (distNorm <= 0);
				if (base3Dnorm) {
					*base3Dnorm = testMe;
					*head3Dnorm = &testMe->norm3D;
				}
                        }
                        else {
                                answer = (distOut <= 0);                        
				if (base3Dnorm) {
					*base3Dnorm = edge->v0;
					*head3Dnorm = &edge->norm3D;
				}
			}
		}
	} while( ++edge < edges + numEdges );
	
	*bestDist2 = bestDistance2;
	return answer;
}

	
