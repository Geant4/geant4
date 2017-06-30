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
//
// $Id: G4GenericTrap.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GenericTrap.cc
//
// Authors:
//   Tatiana Nikitina, CERN; Ivana Hrivnacova, IPN Orsay
//   Adapted from Root Arb8 implementation by Andrei Gheata, CERN 
//
// History:
// 04.08.2011 T.Nikitina - Added SetReferences() and InvertFacets()
//            to CreatePolyhedron() for Visualisation of Boolean
// 03.02.2016 E.Tcherniaev - Revised GetSurfaceArea() and GetCubicVolume(),
//            rewritten GetFaceSurfaceArea(), added GetFaceCubicVolume()      
// 25.09.2016 E.Tcherniaev - Use G4BoundingEnvelope for CalculateExtent(),
//            removed CreateRotatedVertices()
// --------------------------------------------------------------------

#include "G4GenericTrap.hh"

#if !defined(G4GEOM_USE_UGENERICTRAP)

#include <iomanip>

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4TessellatedSolid.hh"
#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "Randomize.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4PolyhedronArbitrary.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

const G4int    G4GenericTrap::fgkNofVertices = 8;
const G4double G4GenericTrap::fgkTolerance = 1E-3;

// --------------------------------------------------------------------

G4GenericTrap::G4GenericTrap( const G4String& name, G4double halfZ,
                              const std::vector<G4TwoVector>&  vertices )
  : G4VSolid(name),
    fRebuildPolyhedron(false),
    fpPolyhedron(0),
    fDz(halfZ),
    fVertices(),
    fIsTwisted(false),
    fTessellatedSolid(0),
    fMinBBoxVector(G4ThreeVector(0,0,0)),
    fMaxBBoxVector(G4ThreeVector(0,0,0)),
    fVisSubdivisions(0),
    fSurfaceArea(0.),
    fCubicVolume(0.)
   
{
  // General constructor
  const G4double min_length=5*1.e-6;
  G4double length = 0.;
  G4int k=0;
  G4String errorDescription = "InvalidSetup in \" ";
  errorDescription += name;
  errorDescription += "\""; 

  halfCarTolerance = kCarTolerance*0.5;

  // Check vertices size

  if ( G4int(vertices.size()) != fgkNofVertices )
  {
    G4Exception("G4GenericTrap::G4GenericTrap()", "GeomSolids0002",
                FatalErrorInArgument, "Number of vertices != 8");
  }            
  
  // Check dZ
  // 
  if (halfZ < kCarTolerance)
  {
     G4Exception("G4GenericTrap::G4GenericTrap()", "GeomSolids0002",
                FatalErrorInArgument, "dZ is too small or negative");
  }           
 
  // Check Ordering and Copy vertices 
  //
  if(CheckOrder(vertices))
  {
    for (G4int i=0; i<fgkNofVertices; ++i) {fVertices.push_back(vertices[i]);}
  }
  else
  { 
    for (G4int i=0; i <4; ++i) {fVertices.push_back(vertices[3-i]);}
    for (G4int i=0; i <4; ++i) {fVertices.push_back(vertices[7-i]);}
  }

   // Check length of segments and Adjust
  // 
  for (G4int j=0; j < 2; j++)
  {
    for (G4int i=1; i<4; ++i)
    {
      k = j*4+i;
      length = (fVertices[k]-fVertices[k-1]).mag();
      if ( ( length < min_length) && ( length > kCarTolerance ) )
      {
        std::ostringstream message;
        message << "Length segment is too small." << G4endl
                << "Distance between " << fVertices[k-1] << " and "
                << fVertices[k] << " is only " << length << " mm !"; 
        G4Exception("G4GenericTrap::G4GenericTrap()", "GeomSolids1001",
                    JustWarning, message, "Vertices will be collapsed.");
        fVertices[k]=fVertices[k-1];
      }
    }
  }

  // Compute Twist
  //
  for( G4int i=0; i<4; i++) { fTwist[i]=0.; }
  fIsTwisted = ComputeIsTwisted();

  // Compute Bounding Box 
  //
  ComputeBBox();
  
  // If not twisted - create tessellated solid 
  // (an alternative implementation for testing)
  //
#ifdef G4TESS_TEST
   if ( !fIsTwisted )  { fTessellatedSolid = CreateTessellatedSolid(); }
#endif
}

// --------------------------------------------------------------------

G4GenericTrap::G4GenericTrap( __void__& a )
  : G4VSolid(a),
    fRebuildPolyhedron(false),
    fpPolyhedron(0),
    halfCarTolerance(0.),
    fDz(0.),
    fVertices(),
    fIsTwisted(false),
    fTessellatedSolid(0),
    fMinBBoxVector(G4ThreeVector(0,0,0)),
    fMaxBBoxVector(G4ThreeVector(0,0,0)),
    fVisSubdivisions(0),
    fSurfaceArea(0.),
    fCubicVolume(0.)
{
  // Fake default constructor - sets only member data and allocates memory
  //                            for usage restricted to object persistency.
}

// --------------------------------------------------------------------

G4GenericTrap::~G4GenericTrap()
{
  // Destructor
  delete fTessellatedSolid;
}

// --------------------------------------------------------------------

G4GenericTrap::G4GenericTrap(const G4GenericTrap& rhs)
  : G4VSolid(rhs),
    fRebuildPolyhedron(false), fpPolyhedron(0),
    halfCarTolerance(rhs.halfCarTolerance),
    fDz(rhs.fDz), fVertices(rhs.fVertices),
    fIsTwisted(rhs.fIsTwisted), fTessellatedSolid(0),
    fMinBBoxVector(rhs.fMinBBoxVector), fMaxBBoxVector(rhs.fMaxBBoxVector),
    fVisSubdivisions(rhs.fVisSubdivisions),
    fSurfaceArea(rhs.fSurfaceArea), fCubicVolume(rhs.fCubicVolume) 
{
   for (size_t i=0; i<4; ++i)  { fTwist[i] = rhs.fTwist[i]; }
#ifdef G4TESS_TEST
   if (rhs.fTessellatedSolid && !fIsTwisted )
   { fTessellatedSolid = CreateTessellatedSolid(); } 
#endif
}

// --------------------------------------------------------------------

G4GenericTrap& G4GenericTrap::operator = (const G4GenericTrap& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   halfCarTolerance = rhs.halfCarTolerance;
   fDz = rhs.fDz; fVertices = rhs.fVertices;
   fIsTwisted = rhs.fIsTwisted; fTessellatedSolid = 0;
   fMinBBoxVector = rhs.fMinBBoxVector; fMaxBBoxVector = rhs.fMaxBBoxVector;
   fVisSubdivisions = rhs.fVisSubdivisions;
   fSurfaceArea = rhs.fSurfaceArea; fCubicVolume = rhs.fCubicVolume;

   for (size_t i=0; i<4; ++i)  { fTwist[i] = rhs.fTwist[i]; }
#ifdef G4TESS_TEST
   if (rhs.fTessellatedSolid && !fIsTwisted )
   { delete fTessellatedSolid; fTessellatedSolid = CreateTessellatedSolid(); } 
#endif
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = 0;

   return *this;
}

// --------------------------------------------------------------------

EInside
G4GenericTrap::InsidePolygone(const G4ThreeVector& p,
                              const std::vector<G4TwoVector>& poly) const 
{
  EInside  in = kInside;
  G4double cross, len2;
  G4int count=0;

  for (G4int i = 0; i < 4; i++)
  {
    G4int j = (i+1) % 4;

    cross = (p.x()-poly[i].x())*(poly[j].y()-poly[i].y())-
            (p.y()-poly[i].y())*(poly[j].x()-poly[i].x());

    len2=(poly[i]-poly[j]).mag2();
    if (len2 > kCarTolerance)
    {
      if(cross*cross<=len2*halfCarTolerance*halfCarTolerance)  // Surface check
      {
        G4double test;

        // Check if p lies between the two extremes of the segment
        //
        G4int iMax;
        G4int iMin;

        if (poly[j].x() > poly[i].x())
        {
          iMax = j;
          iMin = i;
        }
        else {
          iMax = i;
          iMin = j;
        }
        if ( p.x() > poly[iMax].x()+halfCarTolerance
          || p.x() < poly[iMin].x()-halfCarTolerance )
        {
          return kOutside;
        }

        if (poly[j].y() > poly[i].y())
        {
          iMax = j;
          iMin = i;
        }
        else
        {
          iMax = i;
          iMin = j;
        }
        if ( p.y() > poly[iMax].y()+halfCarTolerance
          || p.y() < poly[iMin].y()-halfCarTolerance )
        {
          return kOutside;
        }

        if ( poly[iMax].x() != poly[iMin].x() )
        {
          test = (p.x()-poly[iMin].x())/(poly[iMax].x()-poly[iMin].x())
               * (poly[iMax].y()-poly[iMin].y())+poly[iMin].y();
        }
        else
        {
          test = p.y();
        }

        // Check if point is Inside Segment
        // 
        if( (test>=(poly[iMin].y()-halfCarTolerance))
         && (test<=(poly[iMax].y()+halfCarTolerance)) )
        { 
          return kSurface;
        }
        else
        {
          return kOutside;
        }
      }
      else if (cross<0.)  { return kOutside; }
    }
    else
    {
      count++;
    }
  }

  // All collapsed vertices, Tet like
  //
  if(count==4)
  { 
    if ( (std::fabs(p.x()-poly[0].x())+std::fabs(p.y()-poly[0].y())) > halfCarTolerance )
    {
      in=kOutside;
    }
  }
  return in;
}

// --------------------------------------------------------------------

EInside G4GenericTrap::Inside(const G4ThreeVector& p) const
{
  // Test if point is inside this shape

#ifdef G4TESS_TEST
   if ( fTessellatedSolid )
   { 
     return fTessellatedSolid->Inside(p);
   }
#endif  

  EInside innew=kOutside;
  std::vector<G4TwoVector> xy;
 
  if (std::fabs(p.z()) <= fDz+halfCarTolerance)  // First check Z range
  {
    // Compute intersection between Z plane containing point and the shape
    //
    G4double cf = 0.5*(fDz-p.z())/fDz;
    for (G4int i=0; i<4; i++)
    {
      xy.push_back(fVertices[i+4]+cf*( fVertices[i]-fVertices[i+4]));
    }

    innew=InsidePolygone(p,xy);

    if( (innew==kInside) || (innew==kSurface) )
    { 
      if(std::fabs(p.z()) > fDz-halfCarTolerance)  { innew=kSurface; }
    }
  }
  return innew;    
} 

// --------------------------------------------------------------------

G4ThreeVector G4GenericTrap::SurfaceNormal( const G4ThreeVector& p ) const
{
  // Calculate side nearest to p, and return normal
  // If two sides are equidistant, sum of the Normal is returned

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  {
    return fTessellatedSolid->SurfaceNormal(p);
  }  
#endif   
 
  G4ThreeVector lnorm, sumnorm(0.,0.,0.), apprnorm(0.,0.,1.),
                p0, p1, p2, r1, r2, r3, r4;
  G4int noSurfaces = 0; 
  G4double distxy,distz;
  G4bool zPlusSide=false;

  distz = fDz-std::fabs(p.z());
  if (distz < halfCarTolerance)
  {
    if(p.z()>0)
    {
      zPlusSide=true;
      sumnorm=G4ThreeVector(0,0,1);
    }
    else
    {
      sumnorm=G4ThreeVector(0,0,-1);
    }
    noSurfaces ++;
  } 

  // Check lateral planes
  //
  std:: vector<G4TwoVector> vertices;  
  G4double cf = 0.5*(fDz-p.z())/fDz;
  for (G4int i=0; i<4; i++)
  {
    vertices.push_back(fVertices[i+4]+cf*(fVertices[i]-fVertices[i+4]));
  }

  // Compute distance for lateral planes
  //
  for (G4int q=0; q<4; q++)
  {
    p0=G4ThreeVector(vertices[q].x(),vertices[q].y(),p.z());
    if(zPlusSide)
    {
      p1=G4ThreeVector(fVertices[q].x(),fVertices[q].y(),-fDz);
    }
    else
    {
      p1=G4ThreeVector(fVertices[q+4].x(),fVertices[q+4].y(),fDz); 
    }
    p2=G4ThreeVector(vertices[(q+1)%4].x(),vertices[(q+1)%4].y(),p.z());

    // Collapsed vertices
    //
    if ( (p2-p0).mag2() < kCarTolerance )
    {
      if ( std::fabs(p.z()+fDz) > kCarTolerance )
      {
        p2=G4ThreeVector(fVertices[(q+1)%4].x(),fVertices[(q+1)%4].y(),-fDz);
      }
      else
      {
        p2=G4ThreeVector(fVertices[(q+1)%4+4].x(),fVertices[(q+1)%4+4].y(),fDz);
      }
    }
    lnorm = (p1-p0).cross(p2-p0);
    lnorm = lnorm.unit();
    if(zPlusSide)  { lnorm=-lnorm; }

    // Adjust Normal for Twisted Surface
    //
    if ( (fIsTwisted) && (GetTwistAngle(q)!=0) )
    {
      G4double normP=(p2-p0).mag();
      if(normP)
      {
        G4double proj=(p-p0).dot(p2-p0)/normP;
        if(proj<0)     { proj=0; }
        if(proj>normP) { proj=normP; }
        G4int j=(q+1)%4;
        r1=G4ThreeVector(fVertices[q+4].x(),fVertices[q+4].y(),fDz);
        r2=G4ThreeVector(fVertices[j+4].x(),fVertices[j+4].y(),fDz);
        r3=G4ThreeVector(fVertices[q].x(),fVertices[q].y(),-fDz);
        r4=G4ThreeVector(fVertices[j].x(),fVertices[j].y(),-fDz);
        r1=r1+proj*(r2-r1)/normP;
        r3=r3+proj*(r4-r3)/normP;
        r2=r1-r3;
        r4=r2.cross(p2-p0); r4=r4.unit();
        lnorm=r4;
      }
    }   // End if fIsTwisted

    distxy=std::fabs((p0-p).dot(lnorm));
    if ( distxy<halfCarTolerance )
    {
      noSurfaces ++;

      // Negative sign for Normal is taken for Outside Normal
      //
      sumnorm=sumnorm+lnorm;
    }

    // For ApproxSurfaceNormal
    //
    if (distxy<distz)
    {
      distz=distxy;
      apprnorm=lnorm;
    }      
  }  // End for loop

  // Calculate final Normal, add Normal in the Corners and Touching Sides
  //
  if ( noSurfaces == 0 )
  {
#ifdef G4SPECSDEBUG
    G4Exception("G4GenericTrap::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, "Point p is not on surface !?" );
#endif
    sumnorm=apprnorm;
    // Add Approximative Surface Normal Calculation?
  }
  else if ( noSurfaces == 1 )  { ; }
  else                         { sumnorm = sumnorm.unit(); }

  return sumnorm ; 
}    

// --------------------------------------------------------------------

G4ThreeVector G4GenericTrap::NormalToPlane( const G4ThreeVector& p,
                                            const G4int ipl ) const
{
  // Return normal to given lateral plane ipl

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->SurfaceNormal(p);
  }  
#endif   

  G4ThreeVector lnorm, norm(0.,0.,0.), p0,p1,p2;
 
  G4double  distz = fDz-p.z();
  G4int i=ipl;  // current plane index
 
  G4TwoVector u,v;  
  G4ThreeVector r1,r2,r3,r4;
  G4double cf = 0.5*(fDz-p.z())/fDz;
  G4int j=(i+1)%4;

  u=fVertices[i+4]+cf*(fVertices[i]-fVertices[i+4]);
  v=fVertices[j+4]+cf*(fVertices[j]-fVertices[j+4]);

  // Compute cross product
  //
  p0=G4ThreeVector(u.x(),u.y(),p.z());
      
  if (std::fabs(distz)<halfCarTolerance)
  {
    p1=G4ThreeVector(fVertices[i].x(),fVertices[i].y(),-fDz);
    distz=-1;
  }
  else
  {
    p1=G4ThreeVector(fVertices[i+4].x(),fVertices[i+4].y(),fDz);
  }  
  p2=G4ThreeVector(v.x(),v.y(),p.z());

  // Collapsed vertices
  //
  if ( (p2-p0).mag2() < kCarTolerance )
  {
    if ( std::fabs(p.z()+fDz) > halfCarTolerance )
    {
      p2=G4ThreeVector(fVertices[j].x(),fVertices[j].y(),-fDz);
    }
    else
    {
      p2=G4ThreeVector(fVertices[j+4].x(),fVertices[j+4].y(),fDz);
    }
  }
  lnorm=-(p1-p0).cross(p2-p0);
  if (distz>-halfCarTolerance)  { lnorm=-lnorm.unit(); }
  else                          { lnorm=lnorm.unit();  }
 
  // Adjust Normal for Twisted Surface
  //
  if( (fIsTwisted) && (GetTwistAngle(ipl)!=0) )
  {
    G4double normP=(p2-p0).mag();
    if(normP)
    {
      G4double proj=(p-p0).dot(p2-p0)/normP;
      if (proj<0)     { proj=0; }
      if (proj>normP) { proj=normP; }

      r1=G4ThreeVector(fVertices[i+4].x(),fVertices[i+4].y(),fDz);
      r2=G4ThreeVector(fVertices[j+4].x(),fVertices[j+4].y(),fDz);
      r3=G4ThreeVector(fVertices[i].x(),fVertices[i].y(),-fDz);
      r4=G4ThreeVector(fVertices[j].x(),fVertices[j].y(),-fDz);
      r1=r1+proj*(r2-r1)/normP;
      r3=r3+proj*(r4-r3)/normP;
      r2=r1-r3;
      r4=r2.cross(p2-p0);r4=r4.unit();
      lnorm=r4;
    }
  }  // End if fIsTwisted

  return lnorm;
}    

// --------------------------------------------------------------------

G4double G4GenericTrap::DistToPlane(const G4ThreeVector& p,
                                    const G4ThreeVector& v,
                                    const G4int ipl) const 
{
  // Computes distance to plane ipl :
  // ipl=0 : points 0,4,1,5
  // ipl=1 : points 1,5,2,6
  // ipl=2 : points 2,6,3,7
  // ipl=3 : points 3,7,0,4

  G4double xa,xb,xc,xd,ya,yb,yc,yd;
  
  G4int j = (ipl+1)%4;
  
  xa=fVertices[ipl].x();
  ya=fVertices[ipl].y();
  xb=fVertices[ipl+4].x();
  yb=fVertices[ipl+4].y();
  xc=fVertices[j].x();
  yc=fVertices[j].y();
  xd=fVertices[4+j].x();
  yd=fVertices[4+j].y();
 
  G4double dz2 =0.5/fDz;
  G4double tx1 =dz2*(xb-xa);
  G4double ty1 =dz2*(yb-ya);
  G4double tx2 =dz2*(xd-xc);
  G4double ty2 =dz2*(yd-yc);
  G4double dzp =fDz+p.z();
  G4double xs1 =xa+tx1*dzp;
  G4double ys1 =ya+ty1*dzp;
  G4double xs2 =xc+tx2*dzp;
  G4double ys2 =yc+ty2*dzp;
  G4double dxs =xs2-xs1;
  G4double dys =ys2-ys1;
  G4double dtx =tx2-tx1;
  G4double dty =ty2-ty1;

  G4double a = (dtx*v.y()-dty*v.x()+(tx1*ty2-tx2*ty1)*v.z())*v.z();
  G4double b = dxs*v.y()-dys*v.x()+(dtx*p.y()-dty*p.x()+ty2*xs1-ty1*xs2
             + tx1*ys2-tx2*ys1)*v.z();
  G4double c=dxs*p.y()-dys*p.x()+xs1*ys2-xs2*ys1;
  G4double q=kInfinity;
  G4double x1,x2,y1,y2,xp,yp,zi;

  if (std::fabs(a)<kCarTolerance)
  {           
    if (std::fabs(b)<kCarTolerance)  { return kInfinity; }
    q=-c/b;

    // Check if Point is on the Surface

    if (q>-halfCarTolerance)
    {
      if (q<halfCarTolerance)
      {
        if (NormalToPlane(p,ipl).dot(v)<=0)
          { if(Inside(p) != kOutside) { return 0.; } }
        else
          { return kInfinity; }
      }

      // Check the Intersection
      //
      zi=p.z()+q*v.z();
      if (std::fabs(zi)<fDz)
      {
        x1=xs1+tx1*v.z()*q;
        x2=xs2+tx2*v.z()*q;
        xp=p.x()+q*v.x();
        y1=ys1+ty1*v.z()*q;
        y2=ys2+ty2*v.z()*q;
        yp=p.y()+q*v.y();
        zi = (xp-x1)*(xp-x2)+(yp-y1)*(yp-y2);
        if (zi<=halfCarTolerance)  { return q; }
      }
    }
    return kInfinity;
  }      
  G4double d=b*b-4*a*c;
  if (d>=0)
  {
    if (a>0) { q=0.5*(-b-std::sqrt(d))/a; }
    else     { q=0.5*(-b+std::sqrt(d))/a; }

    // Check if Point is on the Surface
    //
    if (q>-halfCarTolerance)
    {
      if(q<halfCarTolerance)
      {
        if (NormalToPlane(p,ipl).dot(v)<=0)
        {
          if(Inside(p)!= kOutside) { return 0.; }
        }
        else  // Check second root; return kInfinity
        {
          if (a>0) { q=0.5*(-b+std::sqrt(d))/a; }
          else     { q=0.5*(-b-std::sqrt(d))/a; }
          if (q<=halfCarTolerance) { return kInfinity; }
        }
      }
      // Check the Intersection
      //
      zi=p.z()+q*v.z();
      if (std::fabs(zi)<fDz)
      {
        x1=xs1+tx1*v.z()*q;
        x2=xs2+tx2*v.z()*q;
        xp=p.x()+q*v.x();
        y1=ys1+ty1*v.z()*q;
        y2=ys2+ty2*v.z()*q;
        yp=p.y()+q*v.y();
        zi = (xp-x1)*(xp-x2)+(yp-y1)*(yp-y2);
        if (zi<=halfCarTolerance)  { return q; }
      }
    }
    if (a>0)  { q=0.5*(-b+std::sqrt(d))/a; }
    else      { q=0.5*(-b-std::sqrt(d))/a; }

    // Check if Point is on the Surface
    //
    if (q>-halfCarTolerance)
    {
      if(q<halfCarTolerance)
      {
        if (NormalToPlane(p,ipl).dot(v)<=0)
        {
          if(Inside(p) != kOutside)  { return 0.; }
        }
        else   // Check second root; return kInfinity.
        {
          if (a>0) { q=0.5*(-b-std::sqrt(d))/a; }
          else     { q=0.5*(-b+std::sqrt(d))/a; }
          if (q<=halfCarTolerance)  { return kInfinity; }
        }
      }
      // Check the Intersection
      //
      zi=p.z()+q*v.z();
      if (std::fabs(zi)<fDz)
      {
        x1=xs1+tx1*v.z()*q;
        x2=xs2+tx2*v.z()*q;
        xp=p.x()+q*v.x();
        y1=ys1+ty1*v.z()*q;
        y2=ys2+ty2*v.z()*q;
        yp=p.y()+q*v.y();
        zi = (xp-x1)*(xp-x2)+(yp-y1)*(yp-y2);
        if (zi<=halfCarTolerance)  { return q; }
      }
    }
  }
  return kInfinity;
}      

// --------------------------------------------------------------------

G4double G4GenericTrap::DistanceToIn(const G4ThreeVector& p,
                                     const G4ThreeVector& v) const
{
#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->DistanceToIn(p, v);
  }  
#endif
    
  G4double dist[5];
  G4ThreeVector n;

  // Check lateral faces
  //
  G4int i;
  for (i=0; i<4; i++)
  {
    dist[i]=DistToPlane(p, v, i);  
  }

  // Check Z planes
  //
  dist[4]=kInfinity;
  if (std::fabs(p.z())>fDz-halfCarTolerance)
  {
    if (v.z())
    {
      G4ThreeVector pt;
      if (p.z()>0)
      {
        dist[4] = (fDz-p.z())/v.z();
      }
      else
      {   
        dist[4] = (-fDz-p.z())/v.z();
      }   
      if (dist[4]<-halfCarTolerance)
      {
        dist[4]=kInfinity;
      }
      else
      {
        if(dist[4]<halfCarTolerance)
        {
          if(p.z()>0)  { n=G4ThreeVector(0,0,1); }
          else         { n=G4ThreeVector(0,0,-1); }
          if (n.dot(v)<0) { dist[4]=0.; }
          else            { dist[4]=kInfinity; }
        }
        pt=p+dist[4]*v;
        if (Inside(pt)==kOutside)  { dist[4]=kInfinity; }
      }
    }
  }   
  G4double distmin = dist[0];
  for (i=1;i<5;i++)
  {
    if (dist[i] < distmin)  { distmin = dist[i]; }
  }

  if (distmin<halfCarTolerance)  { distmin=0.; }

  return distmin;
}    
  
// --------------------------------------------------------------------

G4double G4GenericTrap::DistanceToIn(const G4ThreeVector& p) const
{
  // Computes the closest distance from given point to this shape

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  {
    return fTessellatedSolid->DistanceToIn(p);
  }  
#endif
 
  G4double safz = std::fabs(p.z())-fDz;
  if(safz<0) { safz=0; }

  G4int iseg;
  G4double safe  = safz;
  G4double safxy = safz;
 
  for (iseg=0; iseg<4; iseg++)
  { 
    safxy = SafetyToFace(p,iseg);
    if (safxy>safe)  { safe=safxy; }
  }

  return safe;
}

// --------------------------------------------------------------------

G4double
G4GenericTrap::SafetyToFace(const G4ThreeVector& p, const G4int iseg) const
{
  // Estimate distance to lateral plane defined by segment iseg in range [0,3]
  // Might be negative: plane seen only from inside

  G4ThreeVector p1,norm;
  G4double safe;

  p1=G4ThreeVector(fVertices[iseg].x(),fVertices[iseg].y(),-fDz);
  
  norm=NormalToPlane(p,iseg);
  safe = (p-p1).dot(norm); // Can be negative
 
  return safe;
}

// --------------------------------------------------------------------

G4double
G4GenericTrap::DistToTriangle(const G4ThreeVector& p,
                              const G4ThreeVector& v, const G4int ipl) const
{
  G4double xa=fVertices[ipl].x();
  G4double ya=fVertices[ipl].y();
  G4double xb=fVertices[ipl+4].x();
  G4double yb=fVertices[ipl+4].y();
  G4int j=(ipl+1)%4;
  G4double xc=fVertices[j].x();
  G4double yc=fVertices[j].y();
  G4double zab=2*fDz;
  G4double zac=0;
  
  if ( (std::fabs(xa-xc)+std::fabs(ya-yc)) < halfCarTolerance )
  {
    xc=fVertices[j+4].x();
    yc=fVertices[j+4].y();
    zac=2*fDz;
    zab=2*fDz;

    //Line case
    //
    if ( (std::fabs(xb-xc)+std::fabs(yb-yc)) < halfCarTolerance )
    {
      return kInfinity;
    }
  }
  G4double a=(yb-ya)*zac-(yc-ya)*zab;
  G4double b=(xc-xa)*zab-(xb-xa)*zac;
  G4double c=(xb-xa)*(yc-ya)-(xc-xa)*(yb-ya);
  G4double d=-xa*a-ya*b+fDz*c;
  G4double t=a*v.x()+b*v.y()+c*v.z();

  if (t!=0)
  {
    t=-(a*p.x()+b*p.y()+c*p.z()+d)/t;
  }
  if ( (t<halfCarTolerance) && (t>-halfCarTolerance) )
  {
    if (NormalToPlane(p,ipl).dot(v)<kCarTolerance)
    {
      t=kInfinity;
    }
    else
    {
      t=0;
    }
  }
  if (Inside(p+v*t) != kSurface)  { t=kInfinity; }
 
  return t;
} 

// --------------------------------------------------------------------

G4double G4GenericTrap::DistanceToOut(const G4ThreeVector& p,
                                      const G4ThreeVector& v,
                                      const G4bool calcNorm,
                                            G4bool* validNorm,
                                            G4ThreeVector* n) const
{
#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->DistanceToOut(p, v, calcNorm, validNorm, n);
  }
#endif

  G4double distmin;
  G4bool lateral_cross = false;
  ESide side = kUndefined;
 
  if (calcNorm)  { *validNorm=true; } // All normals are valid

  if (v.z() < 0)
  {
    distmin=(-fDz-p.z())/v.z();
    if (calcNorm) { side=kMZ; *n=G4ThreeVector(0,0,-1); }
  }
  else
  {
    if (v.z() > 0)
    { 
      distmin = (fDz-p.z())/v.z(); 
      if (calcNorm) { side=kPZ; *n=G4ThreeVector(0,0,1); } 
    }
    else            { distmin = kInfinity; }
  }      

  G4double dz2 =0.5/fDz;
  G4double xa,xb,xc,xd;
  G4double ya,yb,yc,yd;

  for (G4int ipl=0; ipl<4; ipl++)
  {
    G4int j = (ipl+1)%4;
    xa=fVertices[ipl].x();
    ya=fVertices[ipl].y();
    xb=fVertices[ipl+4].x();
    yb=fVertices[ipl+4].y();
    xc=fVertices[j].x();
    yc=fVertices[j].y();
    xd=fVertices[4+j].x();
    yd=fVertices[4+j].y();

    if ( ((std::fabs(xb-xd)+std::fabs(yb-yd))<halfCarTolerance)
      || ((std::fabs(xa-xc)+std::fabs(ya-yc))<halfCarTolerance) )
    {
      G4double q=DistToTriangle(p,v,ipl) ;
      if ( (q>=0) && (q<distmin) )
      {
        distmin=q;
        lateral_cross=true;
        side=ESide(ipl+1);
      }
      continue;
    }
    G4double tx1 =dz2*(xb-xa);
    G4double ty1 =dz2*(yb-ya);
    G4double tx2 =dz2*(xd-xc);
    G4double ty2 =dz2*(yd-yc);
    G4double dzp =fDz+p.z();
    G4double xs1 =xa+tx1*dzp;
    G4double ys1 =ya+ty1*dzp;
    G4double xs2 =xc+tx2*dzp;
    G4double ys2 =yc+ty2*dzp;
    G4double dxs =xs2-xs1;
    G4double dys =ys2-ys1;
    G4double dtx =tx2-tx1;
    G4double dty =ty2-ty1;
    G4double a = (dtx*v.y()-dty*v.x()+(tx1*ty2-tx2*ty1)*v.z())*v.z();
    G4double b = dxs*v.y()-dys*v.x()+(dtx*p.y()-dty*p.x()+ty2*xs1-ty1*xs2
               + tx1*ys2-tx2*ys1)*v.z();
    G4double c=dxs*p.y()-dys*p.x()+xs1*ys2-xs2*ys1;
    G4double q=kInfinity;

    if (std::fabs(a) < kCarTolerance)
    {           
      if (std::fabs(b) < kCarTolerance) { continue; }
      q=-c/b;
         
      // Check for Point on the Surface
      //
      if ((q > -halfCarTolerance) && (q < distmin))
      {
        if (q < halfCarTolerance)
        {
          if (NormalToPlane(p,ipl).dot(v)<0.)  { continue; }
        }
        distmin =q;
        lateral_cross=true;
        side=ESide(ipl+1);
      }   
      continue;
    }
    G4double d=b*b-4*a*c;
    if (d >= 0.)
    {
      if (a > 0)  { q=0.5*(-b-std::sqrt(d))/a; }
      else        { q=0.5*(-b+std::sqrt(d))/a; }

      // Check for Point on the Surface
      //
      if (q > -halfCarTolerance )
      {
        if (q < distmin)
        {
          if(q < halfCarTolerance)
          {
            if (NormalToPlane(p,ipl).dot(v)<0.)  // Check second root
            {
              if (a > 0)  { q=0.5*(-b+std::sqrt(d))/a; }
              else        { q=0.5*(-b-std::sqrt(d))/a; }
              if (( q > halfCarTolerance) && (q < distmin))
              {
                distmin=q;
                lateral_cross = true;
                side=ESide(ipl+1);
              }
              continue;
            }
          }
          distmin = q;
          lateral_cross = true;
          side=ESide(ipl+1);
        }
      }
      else
      {
        if (a > 0)  { q=0.5*(-b+std::sqrt(d))/a; }
        else        { q=0.5*(-b-std::sqrt(d))/a; }

        // Check for Point on the Surface
        //
        if ((q > -halfCarTolerance) && (q < distmin))
        {
          if (q < halfCarTolerance)
          {
            if (NormalToPlane(p,ipl).dot(v)<0.)  // Check second root
            {
              if (a > 0)  { q=0.5*(-b-std::sqrt(d))/a; }
              else        { q=0.5*(-b+std::sqrt(d))/a; }
              if ( ( q > halfCarTolerance) && (q < distmin) )
              {
                distmin=q;
                lateral_cross = true;
                side=ESide(ipl+1);
              }
              continue;
            }  
          }
          distmin =q;
          lateral_cross = true;
          side=ESide(ipl+1);
        }   
      }
    }
  }
  if (!lateral_cross)  // Make sure that track crosses the top or bottom
  {
    if (distmin >= kInfinity)  { distmin=kCarTolerance; }
    G4ThreeVector pt=p+distmin*v;
    
    // Check if propagated point is in the polygon
    //
    G4int i=0;
    if (v.z()>0.) { i=4; }
    std::vector<G4TwoVector> xy;
    for ( G4int j=0; j<4; j++)  { xy.push_back(fVertices[i+j]); }

    // Check Inside
    //
    if (InsidePolygone(pt,xy)==kOutside)
    { 
      if(calcNorm)
      {
        if (v.z()>0) {side= kPZ; *n = G4ThreeVector(0,0,1);}
        else         { side=kMZ; *n = G4ThreeVector(0,0,-1);}
      } 
      return 0.;
    }
    else
    {
      if(v.z()>0) {side=kPZ;}
      else        {side=kMZ;}
    }
  }

  if (calcNorm)
  {
    G4ThreeVector pt=p+v*distmin;     
    switch (side)
    {
      case kXY0:
        *n=NormalToPlane(pt,0);
        break;
      case kXY1:
        *n=NormalToPlane(pt,1);
        break;
      case kXY2:
        *n=NormalToPlane(pt,2);
        break;
      case kXY3:
        *n=NormalToPlane(pt,3);
        break;
      case kPZ:
        *n=G4ThreeVector(0,0,1);
        break;
      case kMZ:
        *n=G4ThreeVector(0,0,-1);
        break;
      default:
        DumpInfo();
        std::ostringstream message;
        G4int oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid." << G4endl
                << "Position:" << G4endl
                << "  p.x() = "   << p.x()/mm << " mm" << G4endl
                << "  p.y() = "   << p.y()/mm << " mm" << G4endl
                << "  p.z() = "   << p.z()/mm << " mm" << G4endl
                << "Direction:" << G4endl
                << "  v.x() = "   << v.x() << G4endl
                << "  v.y() = "   << v.y() << G4endl
                << "  v.z() = "   << v.z() << G4endl
                << "Proposed distance :" << G4endl
                << "  distmin = " << distmin/mm << " mm";
        message.precision(oldprc);
        G4Exception("G4GenericTrap::DistanceToOut(p,v,..)",
                    "GeomSolids1002", JustWarning, message);
        break;
     }
  }
 
  if (distmin<halfCarTolerance)  { distmin=0.; }
 
  return distmin;
}    

// --------------------------------------------------------------------

G4double G4GenericTrap::DistanceToOut(const G4ThreeVector& p) const
{

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->DistanceToOut(p);
  }  
#endif

  G4double safz = fDz-std::fabs(p.z());
  if (safz<0) { safz = 0; }

  G4double safe  = safz;
  G4double safxy = safz;
 
  for (G4int iseg=0; iseg<4; iseg++)
  { 
    safxy = std::fabs(SafetyToFace(p,iseg));
    if (safxy < safe)  { safe = safxy; }
  }

  return safe;
}    

// --------------------------------------------------------------------

void G4GenericTrap::BoundingLimits(G4ThreeVector& pMin,
                                   G4ThreeVector& pMax) const
{
  pMin = GetMinimumBBox();
  pMax = GetMaximumBBox();

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4GenericTrap::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

// --------------------------------------------------------------------

G4bool
G4GenericTrap::CalculateExtent(const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                                     G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  // To build the bounding envelope with plane faces each side face of
  // the trapezoid is subdivided in triangles. Subdivision is done by
  // duplication of vertices in the bases in a way that the envelope be
  // a convex polyhedron (some faces of the envelope can be degenerate)
  //
  G4double dz = GetZHalfLength();
  G4ThreeVectorList baseA(8), baseB(8);
  for (G4int i=0; i<4; ++i)
  {
    G4TwoVector va = GetVertex(i);
    G4TwoVector vb = GetVertex(i+4);
    baseA[2*i].set(va.x(),va.y(),-dz);
    baseB[2*i].set(vb.x(),vb.y(), dz);
  }
  for (G4int i=0; i<4; ++i)
  {
    G4int k1=2*i, k2=(2*i+2)%8;
    G4double ax = (baseA[k2].x()-baseA[k1].x());
    G4double ay = (baseA[k2].y()-baseA[k1].y());
    G4double bx = (baseB[k2].x()-baseB[k1].x());
    G4double by = (baseB[k2].y()-baseB[k1].y());
    G4double znorm = ax*by - ay*bx;
    baseA[k1+1] = (znorm < 0.0) ? baseA[k2] : baseA[k1];
    baseB[k1+1] = (znorm < 0.0) ? baseB[k1] : baseB[k2];
  }

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}
  
// --------------------------------------------------------------------

G4GeometryType G4GenericTrap::GetEntityType() const
{
  return G4String("G4GenericTrap");
}
  
// --------------------------------------------------------------------

G4VSolid* G4GenericTrap::Clone() const
{
  return new G4GenericTrap(*this);
}

// --------------------------------------------------------------------

std::ostream& G4GenericTrap::StreamInfo(std::ostream& os) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " *** \n"
     << "    =================================================== \n"
     << " Solid geometry type: " << GetEntityType()  << G4endl
     << "   half length Z: " << fDz/mm << " mm \n"
     << "   list of vertices:\n";
     
  for ( G4int i=0; i<fgkNofVertices; ++i )
  {
    os << std::setw(5) << "#" << i 
       << "   vx = " << fVertices[i].x()/mm << " mm" 
       << "   vy = " << fVertices[i].y()/mm << " mm" << G4endl;
  }
  os.precision(oldprc);

  return os;
} 

// --------------------------------------------------------------------

G4ThreeVector G4GenericTrap::GetPointOnSurface() const
{

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->GetPointOnSurface();
  }  
#endif

  G4ThreeVector point;
  G4TwoVector u,v,w;
  G4double rand,area,chose,cf,lambda0,lambda1,alfa,beta,zp;
  G4int ipl,j;
 
  std::vector<G4ThreeVector> vertices;
  for (G4int i=0; i<4;i++)
  {
    vertices.push_back(G4ThreeVector(fVertices[i].x(),fVertices[i].y(),-fDz));
  }
  for (G4int i=4; i<8;i++)
  {
    vertices.push_back(G4ThreeVector(fVertices[i].x(),fVertices[i].y(),fDz));
  }

  // Surface Area of Planes(only estimation for twisted)
  //
  G4double Surface0=GetFaceSurfaceArea(vertices[0],vertices[1],
                                       vertices[2],vertices[3]);//-fDz plane 
  G4double Surface1=GetFaceSurfaceArea(vertices[0],vertices[1],
                                       vertices[5],vertices[4]);// Lat plane
  G4double Surface2=GetFaceSurfaceArea(vertices[3],vertices[0],
                                       vertices[4],vertices[7]);// Lat plane
  G4double Surface3=GetFaceSurfaceArea(vertices[2],vertices[3],
                                       vertices[7],vertices[6]);// Lat plane
  G4double Surface4=GetFaceSurfaceArea(vertices[2],vertices[1],
                                       vertices[5],vertices[6]);// Lat plane
  G4double Surface5=GetFaceSurfaceArea(vertices[4],vertices[5],
                                       vertices[6],vertices[7]);// fDz plane
  rand = G4UniformRand();
  area = Surface0+Surface1+Surface2+Surface3+Surface4+Surface5;
  chose = rand*area;

  if ( ( chose < Surface0)
    || ( chose > (Surface0+Surface1+Surface2+Surface3+Surface4)) )
  {                                        // fDz or -fDz Plane
    ipl = G4int(G4UniformRand()*4);
    j = (ipl+1)%4;
    if(chose < Surface0)
    { 
      zp = -fDz;
      u = fVertices[ipl]; v = fVertices[j];
      w = fVertices[(ipl+3)%4]; 
    }
    else
    {
      zp = fDz;
      u = fVertices[ipl+4]; v = fVertices[j+4];
      w = fVertices[(ipl+3)%4+4]; 
    }
    alfa = G4UniformRand();
    beta = G4UniformRand();
    lambda1=alfa*beta;
    lambda0=alfa-lambda1;
    v = v-u;
    w = w-u;
    v = u+lambda0*v+lambda1*w;
  }
  else                                     // Lateral Plane Twisted or Not
  {
    if (chose < Surface0+Surface1) { ipl=0; }
    else if (chose < Surface0+Surface1+Surface2) { ipl=1; }
    else if (chose < Surface0+Surface1+Surface2+Surface3) { ipl=2; }
    else { ipl=3; }
    j = (ipl+1)%4;
    zp = -fDz+G4UniformRand()*2*fDz;
    cf = 0.5*(fDz-zp)/fDz;
    u = fVertices[ipl+4]+cf*( fVertices[ipl]-fVertices[ipl+4]);
    v = fVertices[j+4]+cf*(fVertices[j]-fVertices[j+4]); 
    v = u+(v-u)*G4UniformRand();
  }
  point=G4ThreeVector(v.x(),v.y(),zp);

  return point;
}

// --------------------------------------------------------------------

G4double G4GenericTrap::GetSurfaceArea()
{
  if (fSurfaceArea == 0.0) {
    if(fIsTwisted) {
      fSurfaceArea = G4VSolid::GetSurfaceArea();
    } else {
      // Set vertices
      G4ThreeVector vertix0(fVertices[0].x(),fVertices[0].y(),-fDz);
      G4ThreeVector vertix1(fVertices[1].x(),fVertices[1].y(),-fDz);
      G4ThreeVector vertix2(fVertices[2].x(),fVertices[2].y(),-fDz);
      G4ThreeVector vertix3(fVertices[3].x(),fVertices[3].y(),-fDz);
      G4ThreeVector vertix4(fVertices[4].x(),fVertices[4].y(), fDz);
      G4ThreeVector vertix5(fVertices[5].x(),fVertices[5].y(), fDz);
      G4ThreeVector vertix6(fVertices[6].x(),fVertices[6].y(), fDz);
      G4ThreeVector vertix7(fVertices[7].x(),fVertices[7].y(), fDz);

      // Find Surface Area
      fSurfaceArea = GetFaceSurfaceArea(vertix0,vertix1,vertix2,vertix3)  // -fDz plane
                   + GetFaceSurfaceArea(vertix1,vertix0,vertix4,vertix5)  //  Lat plane
                   + GetFaceSurfaceArea(vertix2,vertix1,vertix5,vertix6)  //  Lat plane 
                   + GetFaceSurfaceArea(vertix3,vertix2,vertix6,vertix7)  //  Lat plane
                   + GetFaceSurfaceArea(vertix0,vertix3,vertix7,vertix4)  //  Lat plane
                   + GetFaceSurfaceArea(vertix7,vertix6,vertix5,vertix4); // +fDz plane 
    }
  }
  return fSurfaceArea;
}

// --------------------------------------------------------------------

G4double G4GenericTrap::GetCubicVolume()
{
  if (fCubicVolume == 0.0) {
    if(fIsTwisted) {
      fCubicVolume = G4VSolid::GetCubicVolume();
    } else {
      // Set vertices
      G4ThreeVector vertix0(fVertices[0].x(),fVertices[0].y(),-fDz);
      G4ThreeVector vertix1(fVertices[1].x(),fVertices[1].y(),-fDz);
      G4ThreeVector vertix2(fVertices[2].x(),fVertices[2].y(),-fDz);
      G4ThreeVector vertix3(fVertices[3].x(),fVertices[3].y(),-fDz);
      G4ThreeVector vertix4(fVertices[4].x(),fVertices[4].y(), fDz);
      G4ThreeVector vertix5(fVertices[5].x(),fVertices[5].y(), fDz);
      G4ThreeVector vertix6(fVertices[6].x(),fVertices[6].y(), fDz);
      G4ThreeVector vertix7(fVertices[7].x(),fVertices[7].y(), fDz);

      // Find Cubic Volume
      fCubicVolume = GetFaceCubicVolume(vertix0,vertix1,vertix2,vertix3)  // -fDz plane
                   + GetFaceCubicVolume(vertix1,vertix0,vertix4,vertix5)  //  Lat plane
                   + GetFaceCubicVolume(vertix2,vertix1,vertix5,vertix6)  //  Lat plane 
                   + GetFaceCubicVolume(vertix3,vertix2,vertix6,vertix7)  //  Lat plane
                   + GetFaceCubicVolume(vertix0,vertix3,vertix7,vertix4)  //  Lat plane
                   + GetFaceCubicVolume(vertix7,vertix6,vertix5,vertix4); // +fDz plane 
    }
  }
  return fCubicVolume;
}

// --------------------------------------------------------------------

G4double G4GenericTrap::GetFaceSurfaceArea(const G4ThreeVector& p0,
                                           const G4ThreeVector& p1, 
                                           const G4ThreeVector& p2,
                                           const G4ThreeVector& p3) const
{
  // Returns area of the facet 
  return (((p2-p0).cross(p3-p1)).mag()) / 2.;
}

// --------------------------------------------------------------------

G4double G4GenericTrap::GetFaceCubicVolume(const G4ThreeVector& p0,
                                           const G4ThreeVector& p1, 
                                           const G4ThreeVector& p2,
                                           const G4ThreeVector& p3) const
{
  // Returns contribution of the facet to the volume of the solid.
  // Orientation of the facet is important, normal should point to outside. 
  return (((p2-p0).cross(p3-p1)).dot(p0)) / 6.;
}

// --------------------------------------------------------------------

G4bool G4GenericTrap::ComputeIsTwisted() 
{
  // Computes tangents of twist angles (angles between projections on XY plane 
  // of corresponding -dz +dz edges). 

  G4bool twisted = false;
  G4double dx1, dy1, dx2, dy2;
  G4int nv = fgkNofVertices/2;

  for ( G4int i=0; i<4; i++ )
  {
    dx1 = fVertices[(i+1)%nv].x()-fVertices[i].x();
    dy1 = fVertices[(i+1)%nv].y()-fVertices[i].y();
    if ( (dx1 == 0) && (dy1 == 0) ) { continue; }

    dx2 = fVertices[nv+(i+1)%nv].x()-fVertices[nv+i].x();
    dy2 = fVertices[nv+(i+1)%nv].y()-fVertices[nv+i].y();

    if ( dx2 == 0 && dy2 == 0 ) { continue; }
    G4double twist_angle = std::fabs(dy1*dx2 - dx1*dy2);        
    if ( twist_angle < fgkTolerance ) { continue; }
    twisted = true;
    SetTwistAngle(i,twist_angle);

    // Check on big angles, potentially navigation problem

    twist_angle = std::acos( (dx1*dx2 + dy1*dy2)
                           / (std::sqrt(dx1*dx1+dy1*dy1)
                             * std::sqrt(dx2*dx2+dy2*dy2)) );
   
    if ( std::fabs(twist_angle) > 0.5*pi+kCarTolerance )
    {
      std::ostringstream message;
      message << "Twisted Angle is bigger than 90 degrees - " << GetName()
              << G4endl
              << "     Potential problem of malformed Solid !" << G4endl
              << "     TwistANGLE = " << twist_angle
              << "*rad  for lateral plane N= " << i;
      G4Exception("G4GenericTrap::ComputeIsTwisted()", "GeomSolids1002",
                  JustWarning, message);
    }
  }

  return twisted;
}

// --------------------------------------------------------------------

G4bool G4GenericTrap::CheckOrder(const std::vector<G4TwoVector>& vertices) const
{
  // Test if the vertices are in a clockwise order, if not reorder them.
  // Also test if they're well defined without crossing opposite segments

  G4bool clockwise_order=true;
  G4double sum1 = 0.;
  G4double sum2 = 0.;
  G4int j;

  for (G4int i=0; i<4; i++)
  {
    j = (i+1)%4;
    sum1 += vertices[i].x()*vertices[j].y() - vertices[j].x()*vertices[i].y();
    sum2 += vertices[i+4].x()*vertices[j+4].y()
          - vertices[j+4].x()*vertices[i+4].y();
  }
  if (sum1*sum2 < -fgkTolerance)
  {
     std::ostringstream message;
     message << "Lower/upper faces defined with opposite clockwise - "
             << GetName();
     G4Exception("G4GenericTrap::CheckOrder()", "GeomSolids0002",
                FatalException, message);
   }
   
   if ((sum1 > 0.)||(sum2 > 0.))
   {
     std::ostringstream message;
     message << "Vertices must be defined in clockwise XY planes - "
             << GetName();
     G4Exception("G4GenericTrap::CheckOrder()", "GeomSolids1001",
                 JustWarning,message, "Re-ordering...");
     clockwise_order = false;
   }

   // Check for illegal crossings
   //
   G4bool illegal_cross = false;
   illegal_cross = IsSegCrossingZ(vertices[0],vertices[4],
                                  vertices[1],vertices[5]);
     
   if (!illegal_cross)
   {
     illegal_cross = IsSegCrossingZ(vertices[2],vertices[6],
                                    vertices[3],vertices[7]);
   }
   // +/- dZ planes
   if (!illegal_cross)
   {
     illegal_cross = IsSegCrossing(vertices[0],vertices[1],
                                   vertices[2],vertices[3]);
   }
   if (!illegal_cross)
   {
     illegal_cross = IsSegCrossing(vertices[0],vertices[3],
                                   vertices[1],vertices[2]);
   }
   if (!illegal_cross)
   {
     illegal_cross = IsSegCrossing(vertices[4],vertices[5],
                                   vertices[6],vertices[7]);
   }
   if (!illegal_cross)
   {
     illegal_cross = IsSegCrossing(vertices[4],vertices[7],
                                   vertices[5],vertices[6]);
   }

   if (illegal_cross)
   {
      std::ostringstream message;
      message << "Malformed polygone with opposite sides - " << GetName();
      G4Exception("G4GenericTrap::CheckOrderAndSetup()",
                  "GeomSolids0002", FatalException, message);
   }
   return clockwise_order;
}

// --------------------------------------------------------------------

void G4GenericTrap::ReorderVertices(std::vector<G4ThreeVector>& vertices) const
{
  // Reorder the vector of vertices 

  std::vector<G4ThreeVector> oldVertices(vertices);

  for ( G4int i=0; i < G4int(oldVertices.size()); ++i )
  {
    vertices[i] = oldVertices[oldVertices.size()-1-i];
  }  
} 
 
// --------------------------------------------------------------------

G4bool
G4GenericTrap::IsSegCrossing(const G4TwoVector& a, const G4TwoVector& b, 
                             const G4TwoVector& c, const G4TwoVector& d) const
{ 
  // Check if segments [A,B] and [C,D] are crossing

  G4bool stand1 = false;
  G4bool stand2 = false;
  G4double dx1,dx2,xm=0.,ym=0.,a1=0.,a2=0.,b1=0.,b2=0.;
  dx1=(b-a).x();
  dx2=(d-c).x();

  if( std::fabs(dx1) < fgkTolerance )  { stand1 = true; }
  if( std::fabs(dx2) < fgkTolerance )  { stand2 = true; }
  if (!stand1)
  {
    a1 = (b.x()*a.y()-a.x()*b.y())/dx1;
    b1 = (b-a).y()/dx1;
  }
  if (!stand2)
  {
    a2 = (d.x()*c.y()-c.x()*d.y())/dx2;
    b2 = (d-c).y()/dx2;
  }   
  if (stand1 && stand2)
  {
    // Segments parallel and vertical
    //
    if (std::fabs(a.x()-c.x())<fgkTolerance)
    {
       // Check if segments are overlapping
       //
       if ( ((c.y()-a.y())*(c.y()-b.y())<-fgkTolerance)
         || ((d.y()-a.y())*(d.y()-b.y())<-fgkTolerance)
         || ((a.y()-c.y())*(a.y()-d.y())<-fgkTolerance)
         || ((b.y()-c.y())*(b.y()-d.y())<-fgkTolerance) )  { return true; }

       return false;
    }
    // Different x values
    //
    return false;
  }
   
  if (stand1)    // First segment vertical
  {
    xm = a.x();
    ym = a2+b2*xm; 
  }
  else
  {
    if (stand2)  // Second segment vertical
    {
      xm = c.x();
      ym = a1+b1*xm;
    }
    else  // Normal crossing
    {
      if (std::fabs(b1-b2) < fgkTolerance)
      {
        // Parallel segments, are they aligned
        //
        if (std::fabs(c.y()-(a1+b1*c.x())) > fgkTolerance)  { return false; }

        // Aligned segments, are they overlapping
        //
        if ( ((c.x()-a.x())*(c.x()-b.x())<-fgkTolerance)
          || ((d.x()-a.x())*(d.x()-b.x())<-fgkTolerance)
          || ((a.x()-c.x())*(a.x()-d.x())<-fgkTolerance)
          || ((b.x()-c.x())*(b.x()-d.x())<-fgkTolerance) )  { return true; }

        return false;
      }
      xm = (a1-a2)/(b2-b1);
      ym = (a1*b2-a2*b1)/(b2-b1);
    }
  }

  // Check if crossing point is both between A,B and C,D
  //
  G4double check = (xm-a.x())*(xm-b.x())+(ym-a.y())*(ym-b.y());
  if (check > -fgkTolerance)  { return false; }
  check = (xm-c.x())*(xm-d.x())+(ym-c.y())*(ym-d.y());
  if (check > -fgkTolerance)  { return false; }

  return true;
}

// --------------------------------------------------------------------

G4bool
G4GenericTrap::IsSegCrossingZ(const G4TwoVector& a, const G4TwoVector& b, 
                              const G4TwoVector& c, const G4TwoVector& d) const
{ 
  // Check if segments [A,B] and [C,D] are crossing when
  // A and C are on -dZ and B and D are on +dZ

  // Calculate the Intersection point between two lines in 3D
  //
  G4ThreeVector temp1,temp2;
  G4ThreeVector v1,v2,p1,p2,p3,p4,dv;
  G4double q,det;
  p1=G4ThreeVector(a.x(),a.y(),-fDz);
  p2=G4ThreeVector(c.x(),c.y(),-fDz);
  p3=G4ThreeVector(b.x(),b.y(),fDz);
  p4=G4ThreeVector(d.x(),d.y(),fDz);
  v1=p3-p1;
  v2=p4-p2;
  dv=p2-p1;

  // In case of Collapsed Vertices No crossing
  //
  if( (std::fabs(dv.x()) < kCarTolerance )&&
      (std::fabs(dv.y()) < kCarTolerance ) )  { return false; }
    
  if( (std::fabs((p4-p3).x()) < kCarTolerance )&&
      (std::fabs((p4-p3).y()) < kCarTolerance ) )  { return false; }
 
  // First estimate if Intersection is possible( if det is 0)
  //
  det = dv.x()*v1.y()*v2.z()+dv.y()*v1.z()*v2.x()
      - dv.x()*v1.z()*v2.y()-dv.y()*v1.x()*v2.z();
   
  if (std::fabs(det)<kCarTolerance)  //Intersection
  {
    temp1 = v1.cross(v2);
    temp2 = (p2-p1).cross(v2);
    if (temp1.dot(temp2) < 0)  { return false; } // intersection negative
    q = temp1.mag();

    if ( q < kCarTolerance )  { return false; }  // parallel lines
    q = ((dv).cross(v2)).mag()/q;
   
    if(q < 1.-kCarTolerance)  { return true; }
  }
  return false;
}

// --------------------------------------------------------------------

G4VFacet*
G4GenericTrap::MakeDownFacet(const std::vector<G4ThreeVector>& fromVertices, 
                             G4int ind1, G4int ind2, G4int ind3) const 
{
  // Create a triangular facet from the polygon points given by indices
  // forming the down side ( the normal goes in -z)
  // Do not create facet if 2 vertices are the same

  if ( (fromVertices[ind1] == fromVertices[ind2]) ||
       (fromVertices[ind2] == fromVertices[ind3]) ||
       (fromVertices[ind1] == fromVertices[ind3]) )  { return 0; }

  std::vector<G4ThreeVector> vertices;
  vertices.push_back(fromVertices[ind1]);
  vertices.push_back(fromVertices[ind2]);
  vertices.push_back(fromVertices[ind3]);
  
  // first vertex most left
  //
  G4ThreeVector cross=(vertices[1]-vertices[0]).cross(vertices[2]-vertices[1]);

  if ( cross.z() > 0.0 )
  {
    // Should not happen, as vertices should have been reordered at this stage

    std::ostringstream message;
    message << "Vertices in wrong order - " << GetName();
    G4Exception("G4GenericTrap::MakeDownFacet", "GeomSolids0002",
                FatalException, message);
  }
  
  return new G4TriangularFacet(vertices[0], vertices[1], vertices[2], ABSOLUTE);
}

// --------------------------------------------------------------------

G4VFacet*
G4GenericTrap::MakeUpFacet(const std::vector<G4ThreeVector>& fromVertices, 
                           G4int ind1, G4int ind2, G4int ind3) const     
{
  // Create a triangular facet from the polygon points given by indices
  // forming the upper side ( z>0 )

  // Do not create facet if 2 vertices are the same
  //
  if ( (fromVertices[ind1] == fromVertices[ind2]) ||
       (fromVertices[ind2] == fromVertices[ind3]) ||
       (fromVertices[ind1] == fromVertices[ind3]) )  { return 0; }

  std::vector<G4ThreeVector> vertices;
  vertices.push_back(fromVertices[ind1]);
  vertices.push_back(fromVertices[ind2]);
  vertices.push_back(fromVertices[ind3]);
  
  // First vertex most left
  //
  G4ThreeVector cross=(vertices[1]-vertices[0]).cross(vertices[2]-vertices[1]);

  if ( cross.z() < 0.0 )
  {
    // Should not happen, as vertices should have been reordered at this stage

    std::ostringstream message;
    message << "Vertices in wrong order - " << GetName();
    G4Exception("G4GenericTrap::MakeUpFacet", "GeomSolids0002",
                FatalException, message);
  }
  
  return new G4TriangularFacet(vertices[0], vertices[1], vertices[2], ABSOLUTE);
}      

// --------------------------------------------------------------------

G4VFacet*
G4GenericTrap::MakeSideFacet(const G4ThreeVector& downVertex0,
                             const G4ThreeVector& downVertex1,
                             const G4ThreeVector& upVertex1,
                             const G4ThreeVector& upVertex0) const      
{
  // Creates a triangular facet from the polygon points given by indices
  // forming the upper side ( z>0 )

  if ( (downVertex0 == downVertex1) && (upVertex0 == upVertex1) )
  {
    return 0;
  }

  if ( downVertex0 == downVertex1 )
  {
    return new G4TriangularFacet(downVertex0, upVertex1, upVertex0, ABSOLUTE);
  }

  if ( upVertex0 == upVertex1 )
  {
    return new G4TriangularFacet(downVertex0, downVertex1, upVertex0, ABSOLUTE);
  }

  return new G4QuadrangularFacet(downVertex0, downVertex1, 
                                 upVertex1, upVertex0, ABSOLUTE);   
}    

// --------------------------------------------------------------------

G4TessellatedSolid* G4GenericTrap::CreateTessellatedSolid() const
{
  // 3D vertices
  //
  G4int nv = fgkNofVertices/2;
  std::vector<G4ThreeVector> downVertices;
  for ( G4int i=0; i<nv; i++ )
  { 
    downVertices.push_back(G4ThreeVector(fVertices[i].x(),
                                         fVertices[i].y(), -fDz));
  }

  std::vector<G4ThreeVector> upVertices;
  for ( G4int i=nv; i<2*nv; i++ )
  { 
    upVertices.push_back(G4ThreeVector(fVertices[i].x(),
                                       fVertices[i].y(), fDz));
  }
                                         
  // Reorder vertices if they are not ordered anti-clock wise
  //
  G4ThreeVector cross 
    = (downVertices[1]-downVertices[0]).cross(downVertices[2]-downVertices[1]);
   G4ThreeVector cross1 
    = (upVertices[1]-upVertices[0]).cross(upVertices[2]-upVertices[1]);
  if ( (cross.z() > 0.0) || (cross1.z() > 0.0) )
  {
    ReorderVertices(downVertices);
    ReorderVertices(upVertices);
  }
    
  G4TessellatedSolid* tessellatedSolid = new G4TessellatedSolid(GetName());
  
  G4VFacet* facet = 0;
  facet = MakeDownFacet(downVertices, 0, 1, 2);
  if (facet)  { tessellatedSolid->AddFacet( facet ); }
  facet = MakeDownFacet(downVertices, 0, 2, 3);
  if (facet)  { tessellatedSolid->AddFacet( facet ); }
  facet = MakeUpFacet(upVertices, 0, 2, 1);
  if (facet)  { tessellatedSolid->AddFacet( facet ); }
  facet = MakeUpFacet(upVertices, 0, 3, 2);
  if (facet)  { tessellatedSolid->AddFacet( facet ); }

  // The quadrangular sides
  //
  for ( G4int i = 0; i < nv; ++i )
  {
    G4int j = (i+1) % nv;
    facet = MakeSideFacet(downVertices[j], downVertices[i], 
                          upVertices[i], upVertices[j]);

    if ( facet )  { tessellatedSolid->AddFacet( facet ); }
  }

  tessellatedSolid->SetSolidClosed(true);

  return tessellatedSolid;
}  

// --------------------------------------------------------------------

void G4GenericTrap::ComputeBBox() 
{
  // Computes bounding box for a shape.

  G4double minX, maxX, minY, maxY;
  minX = maxX = fVertices[0].x();
  minY = maxY = fVertices[0].y();
   
  for (G4int i=1; i< fgkNofVertices; i++)
  {
    if (minX>fVertices[i].x()) { minX=fVertices[i].x(); }
    if (maxX<fVertices[i].x()) { maxX=fVertices[i].x(); }
    if (minY>fVertices[i].y()) { minY=fVertices[i].y(); }
    if (maxY<fVertices[i].y()) { maxY=fVertices[i].y(); }
  }
  fMinBBoxVector = G4ThreeVector(minX,minY,-fDz);
  fMaxBBoxVector = G4ThreeVector(maxX,maxY, fDz);
}

// --------------------------------------------------------------------

G4Polyhedron* G4GenericTrap::GetPolyhedron () const
{

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->GetPolyhedron();
  }
#endif  
  
  if ( (!fpPolyhedron)
    || fRebuildPolyhedron
    || (fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
        fpPolyhedron->GetNumberOfRotationSteps()) )
  {
    G4AutoLock l(&polyhedronMutex);
    delete fpPolyhedron;
    fpPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fpPolyhedron;
}    

// --------------------------------------------------------------------

void G4GenericTrap::DescribeYourselfTo(G4VGraphicsScene& scene) const
{

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->DescribeYourselfTo(scene);
  }
#endif  
  
  scene.AddSolid(*this);
}

// --------------------------------------------------------------------

G4VisExtent G4GenericTrap::GetExtent() const
{
  // Computes bounding vectors for the shape

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->GetExtent();
  } 
#endif
   
  G4ThreeVector minVec = GetMinimumBBox();
  G4ThreeVector maxVec = GetMaximumBBox();
  return G4VisExtent (minVec.x(), maxVec.x(),
                      minVec.y(), maxVec.y(),
                      minVec.z(), maxVec.z()); 
}    

// --------------------------------------------------------------------

G4Polyhedron* G4GenericTrap::CreatePolyhedron() const
{

#ifdef G4TESS_TEST
  if ( fTessellatedSolid )
  { 
    return fTessellatedSolid->CreatePolyhedron();
  }  
#endif 
  
  // Approximation of Twisted Side
  // Construct extra Points, if Twisted Side
  //
  G4PolyhedronArbitrary* polyhedron;
  size_t nVertices, nFacets;

  G4int subdivisions=0;
  G4int i;
  if(fIsTwisted)
  {
    if ( GetVisSubdivisions()!= 0 )
    {
      subdivisions=GetVisSubdivisions();
    }
    else
    {
      // Estimation of Number of Subdivisions for smooth visualisation
      //
      G4double maxTwist=0.;
      for(i=0; i<4; i++)
      {
        if(GetTwistAngle(i)>maxTwist) { maxTwist=GetTwistAngle(i); }
      }

      // Computes bounding vectors for the shape
      //
      G4double Dx,Dy;
      G4ThreeVector minVec = GetMinimumBBox();
      G4ThreeVector maxVec = GetMaximumBBox();
      Dx = 0.5*(maxVec.x()- minVec.y());
      Dy = 0.5*(maxVec.y()- minVec.y());
      if (Dy > Dx)  { Dx=Dy; }
    
      subdivisions=8*G4int(maxTwist/(Dx*Dx*Dx)*fDz);
      if (subdivisions<4)  { subdivisions=4; }
      if (subdivisions>30) { subdivisions=30; }
    }
  }
  G4int sub4=4*subdivisions;
  nVertices = 8+subdivisions*4;
  nFacets = 6+subdivisions*4;
  G4double cf=1./(subdivisions+1);
  polyhedron = new G4PolyhedronArbitrary (nVertices, nFacets);

  // Add Vertex
  //
  for (i=0;i<4;i++)
  {
    polyhedron->AddVertex(G4ThreeVector(fVertices[i].x(),
                                        fVertices[i].y(),-fDz));
  }
  for( i=0;i<subdivisions;i++)
  {
    for(G4int j=0;j<4;j++)
    {
      G4TwoVector u=fVertices[j]+cf*(i+1)*( fVertices[j+4]-fVertices[j]);
      polyhedron->AddVertex(G4ThreeVector(u.x(),u.y(),-fDz+cf*2*fDz*(i+1)));
    }    
  }
  for (i=4;i<8;i++)
  {
    polyhedron->AddVertex(G4ThreeVector(fVertices[i].x(),
                                        fVertices[i].y(),fDz));
  }

  // Add Facets
  //
  polyhedron->AddFacet(1,4,3,2);  //Z-plane
  for (i=0;i<subdivisions+1;i++)
  {
    G4int is=i*4;
    polyhedron->AddFacet(5+is,8+is,4+is,1+is);
    polyhedron->AddFacet(8+is,7+is,3+is,4+is);
    polyhedron->AddFacet(7+is,6+is,2+is,3+is);
    polyhedron->AddFacet(6+is,5+is,1+is,2+is); 
  }
  polyhedron->AddFacet(5+sub4,6+sub4,7+sub4,8+sub4);  //Z-plane

  polyhedron->SetReferences();
  polyhedron->InvertFacets();

  return (G4Polyhedron*) polyhedron;
}

// --------------------------------------------------------------------

#endif
