//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Implementation of G4PolyhedraSide, the face representing
// one segmented side of a Polyhedra
//
// Author: David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4PolyhedraSide.hh"
#include "G4PhysicalConstants.hh"
#include "G4IntersectingCone.hh"
#include "G4ClippablePolygon.hh"
#include "G4AffineTransform.hh"
#include "G4SolidExtentList.hh"
#include "G4GeometryTolerance.hh"

#include "Randomize.hh"

// This new field helps to use the class G4PhSideManager.
//
G4PhSideManager G4PolyhedraSide::subInstanceManager;

// This macro changes the references to fields that are now encapsulated
// in the class G4PhSideData.
//
#define G4MT_phphix ((subInstanceManager.offset[instanceID]).fPhix)
#define G4MT_phphiy ((subInstanceManager.offset[instanceID]).fPhiy)
#define G4MT_phphiz ((subInstanceManager.offset[instanceID]).fPhiz)
#define G4MT_phphik ((subInstanceManager.offset[instanceID]).fPhik)

// Returns the private data instance manager.
//
const G4PhSideManager& G4PolyhedraSide::GetSubInstanceManager()
{
  return subInstanceManager;
}

// Constructor
//
// Values for r1,z1 and r2,z2 should be specified in clockwise
// order in (r,z).
//
G4PolyhedraSide::G4PolyhedraSide( const G4PolyhedraSideRZ* prevRZ,
                                  const G4PolyhedraSideRZ* tail,
                                  const G4PolyhedraSideRZ* head,
                                  const G4PolyhedraSideRZ* nextRZ,
                                        G4int theNumSide, 
                                        G4double thePhiStart, 
                                        G4double thePhiTotal, 
                                        G4bool thePhiIsOpen,
                                        G4bool isAllBehind )
{

  instanceID = subInstanceManager.CreateSubInstance();

  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  G4MT_phphix = 0.0; G4MT_phphiy = 0.0; G4MT_phphiz = 0.0;
  G4MT_phphik = 0.0;

  //
  // Record values
  //
  r[0] = tail->r; z[0] = tail->z;
  r[1] = head->r; z[1] = head->z;
  
  G4double phiTotal;
  
  //
  // Set phi to our convention
  //
  startPhi = thePhiStart;
  while (startPhi < 0.0)    // Loop checking, 13.08.2015, G.Cosmo
    startPhi += twopi;
  
  phiIsOpen = thePhiIsOpen;
  phiTotal = (phiIsOpen) ? thePhiTotal : twopi;
  
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
  G4ThreeVector a1( r[0]*std::cos(phi), r[0]*std::sin(phi), z[0] ),
          b1( r[1]*std::cos(phi), r[1]*std::sin(phi), z[1] ),
          c1( prevRZ->r*std::cos(phi), prevRZ->r*std::sin(phi), prevRZ->z ),
          d1( nextRZ->r*std::cos(phi), nextRZ->r*std::sin(phi), nextRZ->z ),
          a2, b2, c2, d2;
  G4PolyhedraSideEdge *edge = edges;
          
  G4PolyhedraSideVec *vec = vecs;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    //
    // ...this is where we are going
    //
    phi += deltaPhi;
    a2 = G4ThreeVector( r[0]*std::cos(phi), r[0]*std::sin(phi), z[0] );
    b2 = G4ThreeVector( r[1]*std::cos(phi), r[1]*std::sin(phi), z[1] );
    c2 = G4ThreeVector( prevRZ->r*std::cos(phi), prevRZ->r*std::sin(phi), prevRZ->z );
    d2 = G4ThreeVector( nextRZ->r*std::cos(phi), nextRZ->r*std::sin(phi), nextRZ->z );
    
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
    if (vec==vecs)
    {
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
  if (phiIsOpen)
  {
    edge->corner[0] = a2;
    edge->corner[1] = b2;
  }
  else
  {
    vecs[numSide-1].edges[1] = edges;
  }
  
  //
  // Go back and fill in remaining fields in edges
  //
  vec = vecs;
  G4PolyhedraSideVec *prev = vecs+numSide-1;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    edge = vec->edges[0];    // The edge between prev and vec
    
    //
    // Okay: edge normal is average of normals of adjacent faces
    //
    G4ThreeVector eNorm = vec->normal + prev->normal;
    edge->normal = eNorm.unit();  
    
    //
    // Vertex normal is average of norms of adjacent surfaces (all four)
    // However, vec->edgeNorm is unit vector in some direction
    // as the sum of normals of adjacent PolyhedraSide with vec.
    // The normalization used for this vector should be the same
    // for vec and prev.
    //
    eNorm = vec->edgeNorm[0] + prev->edgeNorm[0];
    edge->cornNorm[0] = eNorm.unit();
  
    eNorm = vec->edgeNorm[1] + prev->edgeNorm[1];
    edge->cornNorm[1] = eNorm.unit();
  } while( prev=vec, ++vec < vecs + numSide );
  
  if (phiIsOpen)
  {
    // G4double rFact = std::cos(0.5*deltaPhi);
    //
    // If phi is open, we need to patch up normals of the
    // first and last edges and their corresponding
    // vertices.
    //
    // We use vectors that are in the plane of the
    // face. This should be safe.
    //
    vec = vecs;
    
    G4ThreeVector normvec = vec->edges[0]->corner[0]
                          - vec->edges[0]->corner[1];
    normvec = normvec.cross(vec->normal);
    if (normvec.dot(vec->surfPhi) > 0) normvec = -normvec;

    vec->edges[0]->normal = normvec.unit();
    
    vec->edges[0]->cornNorm[0] = (vec->edges[0]->corner[0]
                                - vec->center).unit();
    vec->edges[0]->cornNorm[1] = (vec->edges[0]->corner[1]
                                - vec->center).unit();
    
    //
    // Repeat for ending phi
    //
    vec = vecs + numSide - 1;
    
    normvec = vec->edges[1]->corner[0] - vec->edges[1]->corner[1];
    normvec = normvec.cross(vec->normal);
    if (normvec.dot(vec->surfPhi) < 0) normvec = -normvec;

    vec->edges[1]->normal = normvec.unit();
    
    vec->edges[1]->cornNorm[0] = (vec->edges[1]->corner[0]
                                - vec->center).unit();
    vec->edges[1]->cornNorm[1] = (vec->edges[1]->corner[1]
                                - vec->center).unit();
  }
  
  //
  // edgeNorm is the factor one multiplies the distance along vector phi
  // on the surface of one of our sides in order to calculate the distance
  // from the edge. (see routine DistanceAway)
  //
  edgeNorm = 1.0/std::sqrt( 1.0 + lenPhi[1]*lenPhi[1] );
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4PolyhedraSide::G4PolyhedraSide( __void__&)
  : startPhi(0.), deltaPhi(0.), endPhi(0.),
    lenRZ(0.), edgeNorm(0.), kCarTolerance(0.), instanceID(0)
{
  r[0] = r[1] = 0.;
  z[0] = z[1] = 0.;
  lenPhi[0] = lenPhi[1] = 0.;
}


// Destructor
//  
G4PolyhedraSide::~G4PolyhedraSide()
{
  delete cone;
  delete [] vecs;
  delete [] edges;
}

// Copy constructor
//
G4PolyhedraSide::G4PolyhedraSide( const G4PolyhedraSide& source )
  : G4VCSGface()
{
  instanceID = subInstanceManager.CreateSubInstance();

  CopyStuff( source );
}


//
// Assignment operator
//
G4PolyhedraSide& G4PolyhedraSide::operator=( const G4PolyhedraSide& source )
{
  if (this == &source) return *this;
  
  delete cone;
  delete [] vecs;
  delete [] edges;
  
  CopyStuff( source );

  return *this;
}

// CopyStuff
//
void G4PolyhedraSide::CopyStuff( const G4PolyhedraSide& source )
{
  //
  // The simple stuff
  //
  numSide    = source.numSide;
  r[0]    = source.r[0];
  r[1]    = source.r[1];
  z[0]    = source.z[0];
  z[1]    = source.z[1];
  startPhi  = source.startPhi;
  deltaPhi  = source.deltaPhi;
  endPhi    = source.endPhi;
  phiIsOpen = source.phiIsOpen;
  allBehind = source.allBehind;
  
  lenRZ     = source.lenRZ;
  lenPhi[0] = source.lenPhi[0];
  lenPhi[1] = source.lenPhi[1];
  edgeNorm  = source.edgeNorm;

  kCarTolerance = source.kCarTolerance;
  fSurfaceArea = source.fSurfaceArea;

  cone = new G4IntersectingCone( *source.cone );

  //
  // Duplicate edges
  //
  G4int  numEdges = phiIsOpen ? numSide+1 : numSide;
  edges = new G4PolyhedraSideEdge[numEdges];
  
  G4PolyhedraSideEdge *edge = edges,
          *sourceEdge = source.edges;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    *edge = *sourceEdge;
  } while( ++sourceEdge, ++edge < edges + numEdges);

  //
  // Duplicate vecs
  //
  vecs = new G4PolyhedraSideVec[numSide];
  
  G4PolyhedraSideVec *vec = vecs,
         *sourceVec = source.vecs;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    *vec = *sourceVec;
    vec->edges[0] = edges + (sourceVec->edges[0] - source.edges);
    vec->edges[1] = edges + (sourceVec->edges[1] - source.edges);
  } while( ++sourceVec, ++vec < vecs + numSide );
}
  
// Intersect
//
// Decide if a line intersects the face.
//
// Arguments:
//  p    = (in) starting point of line segment
//  v    = (in) direction of line segment (assumed a unit vector)
//  A, B    = (in) 2d transform variables (see note top of file)
//  normSign  = (in) desired sign for dot product with normal (see below)
//  surfTolerance  = (in) minimum distance from the surface
//  vecs    = (in) Vector set array
//  distance  = (out) distance to surface furfilling all requirements
//  distFromSurface = (out) distance from the surface
//  thisNormal  = (out) normal vector of the intersecting surface
//
// Return value:
//  true if an intersection is found. Otherwise, output parameters are
//  undefined.
//
// Notes:
// * normSign: if we are "inside" the shape and only want to find out how far
//   to leave the shape, we only want to consider intersections with surfaces in
//   which the trajectory is leaving the shape. Since the normal vectors to the
//   surface always point outwards from the inside, this means we want the dot
//   product of the trajectory direction v and the normal of the side normals[i]
//   to be positive. Thus, we should specify normSign as +1.0. Otherwise, if
//   we are outside and want to go in, normSign should be set to -1.0.
//   Don't set normSign to zero, or you will get no intersections!
//
// * surfTolerance: see notes on argument "surfTolerance" in routine
//   "IntersectSidePlane".
//   ----HOWEVER---- We should *not* apply this surface tolerance if the
//   starting point is not within phi or z of the surface. Specifically,
//   if the starting point p angle in x/y places it on a separate side from the
//   intersection or if the starting point p is outside the z bounds of the
//   segment, surfTolerance must be ignored or we should *always* accept the
//   intersection! 
//   This is simply because the sides do not have infinite extent.
//      
//
G4bool G4PolyhedraSide::Intersect( const G4ThreeVector& p,
                                   const G4ThreeVector& v,  
                                         G4bool outgoing,
                                         G4double surfTolerance,
                                         G4double& distance,
                                         G4double& distFromSurface,
                                         G4ThreeVector& normal,
                                         G4bool& isAllBehind )
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
  G4PolyhedraSideVec* vec = vecs;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
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
    
    if (r[0] > 1/kInfinity && normSign*qa.cross(qc).dot(v) < 0) return false;
    if (r[1] > 1/kInfinity && normSign*qb.cross(qd).dot(v) > 0) return false;
    
    //
    // We allow the face to be slightly behind the trajectory
    // (surface tolerance) only if the point p is within
    // the vicinity of the face
    //
    if (distFromSurface < 0)
    {
      G4ThreeVector ps = p - vec->center; 
      
      G4double rz = ps.dot(vec->surfRZ);
      if (std::fabs(rz) > lenRZ+surfTolerance) return false; 

      G4double pp = ps.dot(vec->surfPhi);
      if (std::fabs(pp) > lenPhi[0]+lenPhi[1]*rz+surfTolerance) return false;
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

// Distance
//
G4double G4PolyhedraSide::Distance( const G4ThreeVector& p, G4bool outgoing )
{
  G4double normSign = outgoing ? -1 : +1;
  
  //
  // Try the closest phi segment first
  //
  G4int iPhi = ClosestPhiSegment( GetPhi(p) );
  
  G4ThreeVector pdotc = p - vecs[iPhi].center;
  G4double normDist = pdotc.dot(vecs[iPhi].normal);
  
  if (normSign*normDist > -0.5*kCarTolerance)
  {
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

// Inside
//
EInside G4PolyhedraSide::Inside( const G4ThreeVector& p,
                                       G4double tolerance, 
                                       G4double* bestDistance )
{
  //
  // Which phi segment is closest to this point?
  //
  G4int iPhi = ClosestPhiSegment( GetPhi(p) );
  
  G4double norm;
  
  //
  // Get distance to this segment
  //
  *bestDistance = DistanceToOneSide( p, vecs[iPhi], &norm );
  
  //
  // Use distance along normal to decide return value
  //
  if ( (std::fabs(norm) > tolerance) || (*bestDistance > 2.0*tolerance) )
    return (norm < 0) ? kInside : kOutside;
  else
    return kSurface;
}

// Normal
//
G4ThreeVector G4PolyhedraSide::Normal( const G4ThreeVector& p,
                                             G4double* bestDistance )
{
  //
  // Which phi segment is closest to this point?
  //
  G4int iPhi = ClosestPhiSegment( GetPhi(p) );

  //
  // Get distance to this segment
  //
  G4double norm;
  *bestDistance = DistanceToOneSide( p, vecs[iPhi], &norm );

  return vecs[iPhi].normal;
}

// Extent
//
G4double G4PolyhedraSide::Extent( const G4ThreeVector axis )
{
  if (axis.perp2() < DBL_MIN)
  {
    //
    // Special case
    //
    return axis.z() < 0 ? -cone->ZLo() : cone->ZHi();
  }

  G4int iPhi, i1, i2;
  G4double best;
  G4ThreeVector* list[4];
  
  //
  // Which phi segment, if any, does the axis belong to
  //
  iPhi = PhiSegment( GetPhi(axis) );
  
  if (iPhi < 0)
  {
    //
    // No phi segment? Check front edge of first side and
    // last edge of second side
    //
    i1 = 0; i2 = numSide-1;
  }
  else
  {
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
  G4ThreeVector** vec = list;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double answer = (*vec)->dot(axis);
    if (answer > best) best = answer;
  } while( ++vec < list+4 );
  
  return best;
}

// CalculateExtent
//
// See notes in G4VCSGface
//
void G4PolyhedraSide::CalculateExtent( const EAxis axis, 
                                       const G4VoxelLimits& voxelLimit,
                                       const G4AffineTransform& transform,
                                             G4SolidExtentList& extentList )
{
  //
  // Loop over all sides
  //
  G4PolyhedraSideVec *vec = vecs;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    //
    // Fill our polygon with the four corners of
    // this side, after the specified transformation
    //
    G4ClippablePolygon polygon;
    
    polygon.AddVertexInOrder(transform.
                             TransformPoint(vec->edges[0]->corner[0]));
    polygon.AddVertexInOrder(transform.
                             TransformPoint(vec->edges[0]->corner[1]));
    polygon.AddVertexInOrder(transform.
                             TransformPoint(vec->edges[1]->corner[1]));
    polygon.AddVertexInOrder(transform.
                             TransformPoint(vec->edges[1]->corner[0]));
    
    //
    // Get extent
    //  
    if (polygon.PartialClip( voxelLimit, axis ))
    {
      //
      // Get dot product of normal along target axis
      //
      polygon.SetNormal( transform.TransformAxis(vec->normal) );

      extentList.AddSurface( polygon );
    }
  } while( ++vec < vecs+numSide );
  
  return;
}

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
//  p    - (in) Point
//  v    - (in) Direction
//  vec    - (in) Description record of the side plane
//  normSign  - (in) Sign (+/- 1) to apply to normal
//  surfTolerance  - (in) Surface tolerance (generally > 0, see below)
//  distance  - (out) Distance along v to intersection
//  distFromSurface - (out) Distance from surface normal
//
// Notes:
//   surfTolerance  - Used to decide if a point is behind the surface,
//        a point is allow to be -surfTolerance behind the
//        surface (as measured along the normal), but *only*
//        if the point is within the r/z bounds + surfTolerance
//        of the segment.
//
G4bool G4PolyhedraSide::IntersectSidePlane( const G4ThreeVector& p,
                                            const G4ThreeVector& v,
                                            const G4PolyhedraSideVec& vec,
                                                  G4double normSign, 
                                                  G4double surfTolerance,
                                                  G4double& distance,
                                                  G4double& distFromSurface )
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
  
  if (atRZ < 0)
  {
    if (r[0]==0) return true;    // Can't miss!
    
    if (atRZ < -lenRZ*1.2) return false;  // Forget it! Missed by a mile.
    
    G4ThreeVector q = p + v;    
    G4ThreeVector qa = q - vec.edges[0]->corner[0],
                  qb = q - vec.edges[1]->corner[0];
    G4ThreeVector qacb = qa.cross(qb);
    if (normSign*qacb.dot(v) < 0) return false;
    
    if (distFromSurface < 0)
    {
      if (atRZ < -lenRZ-surfTolerance) return false;
    }
  }
  else if (atRZ > 0)
  {
    if (r[1]==0) return true;    // Can't miss!
    
    if (atRZ > lenRZ*1.2) return false;  // Missed by a mile
    
    G4ThreeVector q = p + v;    
    G4ThreeVector qa = q - vec.edges[0]->corner[1],
                  qb = q - vec.edges[1]->corner[1];
    G4ThreeVector qacb = qa.cross(qb);
    if (normSign*qacb.dot(v) >= 0) return false;
    
    if (distFromSurface < 0)
    {
      if (atRZ > lenRZ+surfTolerance) return false;
    }
  }

  return true;
}

// LineHitsSegments
//
// Calculate which phi segments a line intersects in three dimensions.
// No check is made as to whether the intersections are within the z bounds of
// the segment.
//
G4int G4PolyhedraSide::LineHitsSegments( const G4ThreeVector& p,
                                         const G4ThreeVector& v,
                                               G4int* i1, G4int* i2 )
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
  *i1 = PhiSegment( std::atan2( p.y() + s1*v.y(), p.x() + s1*v.x() ) );
  if (n==1)
  {
    return (*i1 < 0) ? 0 : 1;
  }
  
  //
  // Try second intersection
  //
  *i2 = PhiSegment( std::atan2( p.y() + s2*v.y(), p.x() + s2*v.x() ) );
  if (*i1 == *i2) return 0;
  
  if (*i1 < 0)
  {
    if (*i2 < 0) return 0;
    *i1 = *i2;
    return 1;
  }

  if (*i2 < 0) return 1;
  
  return 2;
}

// ClosestPhiSegment
//
// Decide which phi segment is closest in phi to the point.
// The result is the same as PhiSegment if there is no phi opening.
//
G4int G4PolyhedraSide::ClosestPhiSegment( G4double phi0 )
{
  G4int iPhi = PhiSegment( phi0 );
  if (iPhi >= 0) return iPhi;
  
  //
  // Boogers! The points falls inside the phi segment.
  // Look for the closest point: the start, or  end
  //
  G4double phi = phi0;
  
  while( phi < startPhi )    // Loop checking, 13.08.2015, G.Cosmo
    phi += twopi;
  G4double d1 = phi-endPhi;

  while( phi > startPhi )    // Loop checking, 13.08.2015, G.Cosmo
    phi -= twopi;
  G4double d2 = startPhi-phi;
  
  return (d2 < d1) ? 0 : numSide-1;
}

// PhiSegment
//
// Decide which phi segment an angle belongs to, counting from zero.
// A value of -1 indicates that the phi value is outside the shape
// (only possible if phiTotal < 360 degrees).
//
G4int G4PolyhedraSide::PhiSegment( G4double phi0 )
{
  //
  // How far are we from phiStart? Come up with a positive answer
  // that is less than 2*PI
  //
  G4double phi = phi0 - startPhi;
  while( phi < 0 )    // Loop checking, 13.08.2015, G.Cosmo
    phi += twopi;
  while( phi > twopi )    // Loop checking, 13.08.2015, G.Cosmo
    phi -= twopi;

  //
  // Divide
  //
  G4int answer = (G4int)(phi/deltaPhi);
  
  if (answer >= numSide)
  {
    if (phiIsOpen)
    {
      return -1;  // Looks like we missed
    }
    else
    {
      answer = numSide-1;  // Probably just roundoff
    }
  }
  
  return answer;
}

// GetPhi
//
// Calculate Phi for a given 3-vector (point), if not already cached for the
// same point, in the attempt to avoid consecutive computation of the same
// quantity
//
G4double G4PolyhedraSide::GetPhi( const G4ThreeVector& p )
{
  G4double val=0.;
  G4ThreeVector vphi(G4MT_phphix, G4MT_phphiy, G4MT_phphiz);

  if (vphi != p)
  {
    val = p.phi();
    G4MT_phphix = p.x(); G4MT_phphiy = p.y(); G4MT_phphiz = p.z();
    G4MT_phphik = val;
  }
  else
  {
    val = G4MT_phphik;
  }
  return val;
}

// DistanceToOneSide
//
// Arguments:
//  p   - (in) Point to check
//  vec   - (in) vector set of this side
//  normDist - (out) distance normal to the side or edge, as appropriate, signed
// Return value = total distance from the side
//
G4double G4PolyhedraSide::DistanceToOneSide( const G4ThreeVector& p,
                                             const G4PolyhedraSideVec& vec,
                                                   G4double* normDist )
{
  G4ThreeVector pct = p - vec.center;
  
  //
  // Get normal distance
  //
  *normDist = vec.normal.dot(pct);

  //
  // Add edge penalty
  //
  return DistanceAway( p, vec, normDist );
}

// DistanceAway
//
// Add distance from side edges, if necessary, to total distance,
// and updates normDist appropriate depending on edge normals.
//
G4double G4PolyhedraSide::DistanceAway( const G4ThreeVector& p,
                                        const G4PolyhedraSideVec& vec,
                                              G4double* normDist )
{
  G4double distOut2;
  G4ThreeVector pct = p - vec.center;
  G4double distFaceNorm = *normDist;
  
  //
  // Okay, are we inside bounds?
  //
  G4double pcDotRZ  = pct.dot(vec.surfRZ);
  G4double pcDotPhi = pct.dot(vec.surfPhi);
  
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
  
  if (pcDotRZ < -lenRZ)
  {
    G4double lenPhiZ = lenPhi[0] - lenRZ*lenPhi[1];
    G4double distOutZ = pcDotRZ+lenRZ;
    //
    // Below in RZ
    //
    if (pcDotPhi < -lenPhiZ)
    {
      //
      // ...and below in phi. Find distance to point (A)
      //
      G4double distOutPhi = pcDotPhi+lenPhiZ;
      distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
      G4ThreeVector pa = p - vec.edges[0]->corner[0];
      *normDist = pa.dot(vec.edges[0]->cornNorm[0]);
    }
    else if (pcDotPhi > lenPhiZ)
    {
      //
      // ...and above in phi. Find distance to point (B)
      //
      G4double distOutPhi = pcDotPhi-lenPhiZ;
      distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
      G4ThreeVector pb = p - vec.edges[1]->corner[0];
      *normDist = pb.dot(vec.edges[1]->cornNorm[0]);
    }
    else
    {
      //
      // ...and inside in phi. Find distance to line (C)
      //
      G4ThreeVector pa = p - vec.edges[0]->corner[0];
      distOut2 = distOutZ*distOutZ;
      *normDist = pa.dot(vec.edgeNorm[0]);
    }
  }
  else if (pcDotRZ > lenRZ)
  {
    G4double lenPhiZ = lenPhi[0] + lenRZ*lenPhi[1];
    G4double distOutZ = pcDotRZ-lenRZ;
    //
    // Above in RZ
    //
    if (pcDotPhi < -lenPhiZ)
    {
      //
      // ...and below in phi. Find distance to point (D)
      //
      G4double distOutPhi = pcDotPhi+lenPhiZ;
      distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
      G4ThreeVector pd = p - vec.edges[0]->corner[1];
      *normDist = pd.dot(vec.edges[0]->cornNorm[1]);
    }
    else if (pcDotPhi > lenPhiZ)
    {
      //
      // ...and above in phi. Find distance to point (E)
      //
      G4double distOutPhi = pcDotPhi-lenPhiZ;
      distOut2 = distOutPhi*distOutPhi + distOutZ*distOutZ;
      G4ThreeVector pe = p - vec.edges[1]->corner[1];
      *normDist = pe.dot(vec.edges[1]->cornNorm[1]);
    }
    else
    {
      //
      // ...and inside in phi. Find distance to line (F)
      //
      distOut2 = distOutZ*distOutZ;
      G4ThreeVector pd = p - vec.edges[0]->corner[1];
      *normDist = pd.dot(vec.edgeNorm[1]);
    }
  }
  else
  {
    G4double lenPhiZ = lenPhi[0] + pcDotRZ*lenPhi[1];
    //
    // We are inside RZ bounds
    // 
    if (pcDotPhi < -lenPhiZ)
    {
      //
      // ...and below in phi. Find distance to line (G)
      //
      G4double distOut = edgeNorm*(pcDotPhi+lenPhiZ);
      distOut2 = distOut*distOut;
      G4ThreeVector pd = p - vec.edges[0]->corner[1];
      *normDist = pd.dot(vec.edges[0]->normal);
    }
    else if (pcDotPhi > lenPhiZ)
    {
      //
      // ...and above in phi. Find distance to line (H)
      //
      G4double distOut = edgeNorm*(pcDotPhi-lenPhiZ);
      distOut2 = distOut*distOut;
      G4ThreeVector pe = p - vec.edges[1]->corner[1];
      *normDist = pe.dot(vec.edges[1]->normal);
    }
    else
    {
      //
      // Inside bounds! No penalty.
      //
      return std::fabs(distFaceNorm);
    }
  }
  return std::sqrt( distFaceNorm*distFaceNorm + distOut2 );
}

// Calculation of surface area of a triangle. 
// At the same time a random point in the triangle is given
//
G4double G4PolyhedraSide::SurfaceTriangle( G4ThreeVector p1,
                                           G4ThreeVector p2,
                                           G4ThreeVector p3,
                                           G4ThreeVector* p4 )
{
  G4ThreeVector v, w;
  
  v = p3 - p1;
  w = p1 - p2;
  G4double lambda1 = G4UniformRand();
  G4double lambda2 = lambda1*G4UniformRand();
 
  *p4=p2 + lambda1*w + lambda2*v;
  return 0.5*(v.cross(w)).mag();
}

// GetPointOnPlane
//
// Auxiliary method for GetPointOnSurface()
//
G4ThreeVector
G4PolyhedraSide::GetPointOnPlane( G4ThreeVector p0, G4ThreeVector p1, 
                                  G4ThreeVector p2, G4ThreeVector p3,
                                  G4double* Area )
{
  G4double chose,aOne,aTwo;
  G4ThreeVector point1,point2;
  aOne = SurfaceTriangle(p0,p1,p2,&point1);
  aTwo = SurfaceTriangle(p2,p3,p0,&point2);
  *Area= aOne+aTwo;

  chose = G4UniformRand()*(aOne+aTwo);
  if( (chose>=0.) && (chose < aOne) )
  {
   return (point1);    
  }
  return (point2);
}

// SurfaceArea()
//
G4double G4PolyhedraSide::SurfaceArea()
{
  if( fSurfaceArea==0. )
  { 
    // Define the variables
    //
    G4double area,areas;
    G4ThreeVector point1;
    G4ThreeVector v1,v2,v3,v4; 
    G4PolyhedraSideVec* vec = vecs;
    areas=0.;

    // Do a loop on all SideEdge
    //
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      // Define 4points for a Plane or Triangle
      //
      v1=vec->edges[0]->corner[0];
      v2=vec->edges[0]->corner[1];
      v3=vec->edges[1]->corner[1];
      v4=vec->edges[1]->corner[0];
      point1=GetPointOnPlane(v1,v2,v3,v4,&area);
      areas+=area;
    } while( ++vec < vecs + numSide);

    fSurfaceArea=areas;
  }
  return fSurfaceArea;
}

// GetPointOnFace()
//
G4ThreeVector G4PolyhedraSide::GetPointOnFace()
{
  // Define the variables
  //
  std::vector<G4double>areas;
  std::vector<G4ThreeVector>points;
  G4double area=0.;
  G4double result1;
  G4ThreeVector point1;
  G4ThreeVector v1,v2,v3,v4; 
  G4PolyhedraSideVec* vec = vecs;

  // Do a loop on all SideEdge
  //
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    // Define 4points for a Plane or Triangle
    //
    v1=vec->edges[0]->corner[0];
    v2=vec->edges[0]->corner[1];
    v3=vec->edges[1]->corner[1];
    v4=vec->edges[1]->corner[0];
    point1=GetPointOnPlane(v1,v2,v3,v4,&result1);
    points.push_back(point1);
    areas.push_back(result1);
    area+=result1;
  } while( ++vec < vecs+numSide );

  // Choose randomly one of the surfaces and point on it
  //
  G4double chose = area*G4UniformRand();
  G4double Achose1=0., Achose2=0.;
  G4int i=0;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    Achose2+=areas[i];
    if(chose>=Achose1 && chose<Achose2)
    {
      point1=points[i] ; break;     
    }
    ++i; Achose1=Achose2;
  } while( i<numSide );
 
  return point1;
}
