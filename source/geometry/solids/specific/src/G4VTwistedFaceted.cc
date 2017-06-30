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
// $Id: G4VTwistedFaceted.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4VTwistedFaceted.cc
//
// Author:
//
//   04-Nov-2004 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include "G4VTwistedFaceted.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4SolidExtentList.hh"
#include "G4ClippablePolygon.hh"
#include "G4VPVParameterisation.hh"
#include "G4GeometryTolerance.hh"
#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"

#include "Randomize.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

//=====================================================================
//* constructors ------------------------------------------------------

G4VTwistedFaceted::
G4VTwistedFaceted( const G4String &pname,     // Name of instance
                         G4double  PhiTwist,  // twist angle
                         G4double  pDz,       // half z length
                         G4double  pTheta, // direction between end planes
                         G4double  pPhi,   // defined by polar and azim. angles
                         G4double  pDy1,   // half y length at -pDz
                         G4double  pDx1,   // half x length at -pDz,-pDy
                         G4double  pDx2,   // half x length at -pDz,+pDy
                         G4double  pDy2,   // half y length at +pDz
                         G4double  pDx3,   // half x length at +pDz,-pDy
                         G4double  pDx4,   // half x length at +pDz,+pDy
                         G4double  pAlph   // tilt angle 
                 )
  : G4VSolid(pname), fRebuildPolyhedron(false), fpPolyhedron(0),
    fLowerEndcap(0), fUpperEndcap(0), fSide0(0),
    fSide90(0), fSide180(0), fSide270(0),
    fSurfaceArea(0.)
{

  G4double pDytmp ;
  G4double fDxUp ;
  G4double fDxDown ;

  fDx1 = pDx1 ;
  fDx2 = pDx2 ;
  fDx3 = pDx3 ;
  fDx4 = pDx4 ;
  fDy1 = pDy1 ;
  fDy2 = pDy2 ;
  fDz  = pDz  ;

  G4double kAngTolerance
    = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  // maximum values
  //
  fDxDown = ( fDx1 > fDx2 ? fDx1 : fDx2 ) ;
  fDxUp   = ( fDx3 > fDx4 ? fDx3 : fDx4 ) ;
  fDx     = ( fDxUp > fDxDown ? fDxUp : fDxDown ) ;
  fDy     = ( fDy1 > fDy2 ? fDy1 : fDy2 ) ;
  
  // planarity check
  //
  if ( fDx1 != fDx2 && fDx3 != fDx4 )
  {
    pDytmp = fDy1 * ( fDx3 - fDx4 ) / ( fDx1 - fDx2 ) ;
    if ( std::fabs(pDytmp - fDy2) > kCarTolerance )
    {
      std::ostringstream message;
      message << "Not planar surface in untwisted Trapezoid: "
              << GetName() << G4endl
              << "fDy2 is " << fDy2 << " but should be "
              << pDytmp << ".";
      G4Exception("G4VTwistedFaceted::G4VTwistedFaceted()", "GeomSolids0002",
                  FatalErrorInArgument, message);
    }
  }

#ifdef G4TWISTDEBUG
  if ( fDx1 == fDx2 && fDx3 == fDx4 )
  { 
      G4cout << "Trapezoid is a box" << G4endl ;
  }
  
#endif

  if ( (  fDx1 == fDx2 && fDx3 != fDx4 ) || ( fDx1 != fDx2 && fDx3 == fDx4 ) )
  {
    std::ostringstream message;
    message << "Not planar surface in untwisted Trapezoid: "
            << GetName() << G4endl
            << "One endcap is rectangular, the other is a trapezoid." << G4endl
            << "For planarity reasons they have to be rectangles or trapezoids "
            << "on both sides.";
    G4Exception("G4VTwistedFaceted::G4VTwistedFaceted()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  // twist angle
  //
  fPhiTwist = PhiTwist ;

  // tilt angle
  //
  fAlph = pAlph ; 
  fTAlph = std::tan(fAlph) ;
  
  fTheta = pTheta ;
  fPhi   = pPhi ;

  // dx in surface equation
  //
  fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi)  ;

  // dy in surface equation
  //
  fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi)  ;

  if  ( ! ( ( fDx1  > 2*kCarTolerance)
         && ( fDx2  > 2*kCarTolerance)
         && ( fDx3  > 2*kCarTolerance)
         && ( fDx4  > 2*kCarTolerance)
         && ( fDy1  > 2*kCarTolerance)
         && ( fDy2  > 2*kCarTolerance)
         && ( fDz   > 2*kCarTolerance) 
         && ( std::fabs(fPhiTwist) > 2*kAngTolerance )
         && ( std::fabs(fPhiTwist) < pi/2 )
         && ( std::fabs(fAlph) < pi/2 )
         && ( fTheta < pi/2 && fTheta >= 0 ) )
      )
  {
    std::ostringstream message;
    message << "Invalid dimensions. Too small, or twist angle too big: "
            << GetName() << G4endl
            << "fDx 1-4 = " << fDx1/cm << ", " << fDx2/cm << ", "
            << fDx3/cm << ", " << fDx4/cm << " cm" << G4endl 
            << "fDy 1-2 = " << fDy1/cm << ", " << fDy2/cm << ", "
            << " cm" << G4endl 
            << "fDz = " << fDz/cm << " cm" << G4endl 
            << " twistangle " << fPhiTwist/deg << " deg" << G4endl 
            << " phi,theta = " << fPhi/deg << ", "  << fTheta/deg << " deg";
    G4Exception("G4TwistedTrap::G4VTwistedFaceted()",
                "GeomSolids0002", FatalErrorInArgument, message);
  }
  CreateSurfaces();
  fCubicVolume = 2 * fDz * ( ( fDx1 + fDx2 ) * fDy1 + ( fDx3 + fDx4 ) * fDy2 );
}


//=====================================================================
//* Fake default constructor ------------------------------------------

G4VTwistedFaceted::G4VTwistedFaceted( __void__& a )
  : G4VSolid(a), fRebuildPolyhedron(false), fpPolyhedron(0),
    fTheta(0.), fPhi(0.), fDy1(0.),
    fDx1(0.), fDx2(0.), fDy2(0.), fDx3(0.), fDx4(0.),
    fDz(0.), fDx(0.), fDy(0.), fAlph(0.),
    fTAlph(0.), fdeltaX(0.), fdeltaY(0.), fPhiTwist(0.),
    fLowerEndcap(0), fUpperEndcap(0), fSide0(0), fSide90(0), fSide180(0),
    fSide270(0), fCubicVolume(0.), fSurfaceArea(0.)
{
}


//=====================================================================
//* destructor --------------------------------------------------------

G4VTwistedFaceted::~G4VTwistedFaceted()
{
  if (fLowerEndcap) { delete fLowerEndcap ; }
  if (fUpperEndcap) { delete fUpperEndcap ; }

  if (fSide0)       { delete fSide0   ; }
  if (fSide90)      { delete fSide90  ; }
  if (fSide180)     { delete fSide180 ; }
  if (fSide270)     { delete fSide270 ; }
  if (fpPolyhedron) { delete fpPolyhedron; fpPolyhedron = 0; }
}


//=====================================================================
//* Copy constructor --------------------------------------------------

G4VTwistedFaceted::G4VTwistedFaceted(const G4VTwistedFaceted& rhs)
  : G4VSolid(rhs), fRebuildPolyhedron(false), fpPolyhedron(0),
    fTheta(rhs.fTheta), fPhi(rhs.fPhi),
    fDy1(rhs.fDy1), fDx1(rhs.fDx1), fDx2(rhs.fDx2), fDy2(rhs.fDy2),
    fDx3(rhs.fDx3), fDx4(rhs.fDx4), fDz(rhs.fDz), fDx(rhs.fDx), fDy(rhs.fDy),
    fAlph(rhs.fAlph), fTAlph(rhs.fTAlph), fdeltaX(rhs.fdeltaX),
    fdeltaY(rhs.fdeltaY), fPhiTwist(rhs.fPhiTwist), fLowerEndcap(0),
    fUpperEndcap(0), fSide0(0), fSide90(0), fSide180(0), fSide270(0),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    fLastInside(rhs.fLastInside), fLastNormal(rhs.fLastNormal),
    fLastDistanceToIn(rhs.fLastDistanceToIn),
    fLastDistanceToOut(rhs.fLastDistanceToOut),
    fLastDistanceToInWithV(rhs.fLastDistanceToInWithV),
    fLastDistanceToOutWithV(rhs.fLastDistanceToOutWithV)
{
  CreateSurfaces();
}


//=====================================================================
//* Assignment operator -----------------------------------------------

G4VTwistedFaceted& G4VTwistedFaceted::operator = (const G4VTwistedFaceted& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   fTheta = rhs.fTheta; fPhi = rhs.fPhi;
   fDy1= rhs.fDy1; fDx1= rhs.fDx1; fDx2= rhs.fDx2; fDy2= rhs.fDy2;
   fDx3= rhs.fDx3; fDx4= rhs.fDx4; fDz= rhs.fDz; fDx= rhs.fDx; fDy= rhs.fDy;
   fAlph= rhs.fAlph; fTAlph= rhs.fTAlph; fdeltaX= rhs.fdeltaX;
   fdeltaY= rhs.fdeltaY; fPhiTwist= rhs.fPhiTwist; fLowerEndcap= 0;
   fUpperEndcap= 0; fSide0= 0; fSide90= 0; fSide180= 0; fSide270= 0;
   fCubicVolume= rhs.fCubicVolume; fSurfaceArea= rhs.fSurfaceArea;
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron= 0;
   fLastInside= rhs.fLastInside; fLastNormal= rhs.fLastNormal;
   fLastDistanceToIn= rhs.fLastDistanceToIn;
   fLastDistanceToOut= rhs.fLastDistanceToOut;
   fLastDistanceToInWithV= rhs.fLastDistanceToInWithV;
   fLastDistanceToOutWithV= rhs.fLastDistanceToOutWithV;
 
   CreateSurfaces();

   return *this;
}


//=====================================================================
//* ComputeDimensions -------------------------------------------------

void G4VTwistedFaceted::ComputeDimensions(G4VPVParameterisation* ,
                                          const G4int ,
                                          const G4VPhysicalVolume* )
{
  G4Exception("G4VTwistedFaceted::ComputeDimensions()",
              "GeomSolids0001", FatalException,
              "G4VTwistedFaceted does not support Parameterisation.");
}


//=====================================================================
//* Extent ------------------------------------------------------------

void G4VTwistedFaceted::BoundingLimits(G4ThreeVector &pMin,
                                       G4ThreeVector &pMax) const
{
  G4double maxRad = std::sqrt(fDx*fDx + fDy*fDy);
  pMin.set(-maxRad,-maxRad,-fDz);
  pMax.set( maxRad, maxRad, fDz);
}


//=====================================================================
//* CalculateExtent ---------------------------------------------------

G4bool
G4VTwistedFaceted::CalculateExtent( const EAxis              pAxis,
                                    const G4VoxelLimits     &pVoxelLimit,
                                    const G4AffineTransform &pTransform,
                                          G4double          &pMin,
                                          G4double          &pMax ) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}


//=====================================================================
//* Inside ------------------------------------------------------------

EInside G4VTwistedFaceted::Inside(const G4ThreeVector& p) const
{

   G4ThreeVector *tmpp;
   EInside       *tmpin;
   if (fLastInside.p == p) {
     return fLastInside.inside;
   } else {
      tmpp      = const_cast<G4ThreeVector*>(&(fLastInside.p));
      tmpin     = const_cast<EInside*>(&(fLastInside.inside));
      tmpp->set(p.x(), p.y(), p.z());
   }

   *tmpin = kOutside ;

   G4double phi = p.z()/(2*fDz) * fPhiTwist ;  // rotate the point to z=0
   G4double cphi = std::cos(-phi) ;
   G4double sphi = std::sin(-phi) ;

   G4double px  = p.x() + fdeltaX * ( -phi/fPhiTwist) ;  // shift
   G4double py  = p.y() + fdeltaY * ( -phi/fPhiTwist) ;
   G4double pz  = p.z() ;

   G4double posx = px * cphi - py * sphi   ;  // rotation
   G4double posy = px * sphi + py * cphi   ;
   G4double posz = pz  ;

   G4double xMin = Xcoef(posy,phi,fTAlph) - 2*Xcoef(posy,phi,0.) ; 
   G4double xMax = Xcoef(posy,phi,fTAlph) ;  

   G4double yMax = GetValueB(phi)/2. ;  // b(phi)/2 is limit
   G4double yMin = -yMax ;

#ifdef G4TWISTDEBUG

   G4cout << "inside called: p = " << p << G4endl ; 
   G4cout << "fDx1 = " << fDx1 << G4endl ;
   G4cout << "fDx2 = " << fDx2 << G4endl ;
   G4cout << "fDx3 = " << fDx3 << G4endl ;
   G4cout << "fDx4 = " << fDx4 << G4endl ;

   G4cout << "fDy1 = " << fDy1 << G4endl ;
   G4cout << "fDy2 = " << fDy2 << G4endl ;

   G4cout << "fDz  = " << fDz  << G4endl ;

   G4cout << "Tilt angle alpha = " << fAlph << G4endl ;
   G4cout << "phi,theta = " << fPhi << " , " << fTheta << G4endl ;

   G4cout << "Twist angle = " << fPhiTwist << G4endl ;

   G4cout << "posx = " << posx << G4endl ;
   G4cout << "posy = " << posy << G4endl ;
   G4cout << "xMin = " << xMin << G4endl ;
   G4cout << "xMax = " << xMax << G4endl ;
   G4cout << "yMin = " << yMin << G4endl ;
   G4cout << "yMax = " << yMax << G4endl ;

#endif 


  if ( posx <= xMax - kCarTolerance*0.5
    && posx >= xMin + kCarTolerance*0.5 )
  {
    if ( posy <= yMax - kCarTolerance*0.5
      && posy >= yMin + kCarTolerance*0.5 )
    {
      if      (std::fabs(posz) <= fDz - kCarTolerance*0.5 ) *tmpin = kInside ;
      else if (std::fabs(posz) <= fDz + kCarTolerance*0.5 ) *tmpin = kSurface ;
    }
    else if ( posy <= yMax + kCarTolerance*0.5
           && posy >= yMin - kCarTolerance*0.5 )
    {
      if (std::fabs(posz) <= fDz + kCarTolerance*0.5 ) *tmpin = kSurface ;
    }
  }
  else if ( posx <= xMax + kCarTolerance*0.5
         && posx >= xMin - kCarTolerance*0.5 )
  {
    if ( posy <= yMax + kCarTolerance*0.5
      && posy >= yMin - kCarTolerance*0.5 )
    {
      if (std::fabs(posz) <= fDz + kCarTolerance*0.5) *tmpin = kSurface ;
    }
  }

#ifdef G4TWISTDEBUG
  G4cout << "inside = " << fLastInside.inside << G4endl ;
#endif

  return fLastInside.inside;

}


//=====================================================================
//* SurfaceNormal -----------------------------------------------------

G4ThreeVector G4VTwistedFaceted::SurfaceNormal(const G4ThreeVector& p) const
{
   //
   // return the normal unit vector to the Hyperbolical Surface at a point 
   // p on (or nearly on) the surface
   //
   // Which of the three or four surfaces are we closest to?
   //

   if (fLastNormal.p == p)
   {
     return fLastNormal.vec;
   } 
   
   G4ThreeVector *tmpp       = const_cast<G4ThreeVector*>(&(fLastNormal.p));
   G4ThreeVector *tmpnormal  = const_cast<G4ThreeVector*>(&(fLastNormal.vec));
   G4VTwistSurface    **tmpsurface = const_cast<G4VTwistSurface**>(fLastNormal.surface);
   tmpp->set(p.x(), p.y(), p.z());

   G4double      distance = kInfinity;

   G4VTwistSurface *surfaces[6];

   surfaces[0] = fSide0 ;
   surfaces[1] = fSide90 ;
   surfaces[2] = fSide180 ;
   surfaces[3] = fSide270 ;
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;

   G4ThreeVector xx;
   G4ThreeVector bestxx;
   G4int i;
   G4int besti = -1;
   for (i=0; i< 6; i++)
   {
      G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
      if (tmpdistance < distance)
      {
         distance = tmpdistance;
         bestxx = xx;
         besti = i; 
      }
   }

   tmpsurface[0] = surfaces[besti];
   *tmpnormal = tmpsurface[0]->GetNormal(bestxx, true);
   
   return fLastNormal.vec;
}


//=====================================================================
//* DistanceToIn (p, v) -----------------------------------------------

G4double G4VTwistedFaceted::DistanceToIn (const G4ThreeVector& p,
                                          const G4ThreeVector& v ) const
{

   // DistanceToIn (p, v):
   // Calculate distance to surface of shape from `outside' 
   // along with the v, allowing for tolerance.
   // The function returns kInfinity if no intersection or
   // just grazing within tolerance.

   //
   // checking last value
   //
   
   G4ThreeVector *tmpp;
   G4ThreeVector *tmpv;
   G4double      *tmpdist;
   if (fLastDistanceToInWithV.p == p && fLastDistanceToInWithV.vec == v)
   {
      return fLastDistanceToIn.value;
   }
   else
   {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToInWithV.p));
      tmpv    = const_cast<G4ThreeVector*>(&(fLastDistanceToInWithV.vec));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToInWithV.value));
      tmpp->set(p.x(), p.y(), p.z());
      tmpv->set(v.x(), v.y(), v.z());
   }

   //
   // Calculate DistanceToIn(p,v)
   //
   
   EInside currentside = Inside(p);

   if (currentside == kInside)
   {
   }
   else if (currentside == kSurface)
   {
     // particle is just on a boundary.
     // if the particle is entering to the volume, return 0
     //
     G4ThreeVector normal = SurfaceNormal(p);
     if (normal*v < 0)
     {
       *tmpdist = 0;
       return fLastDistanceToInWithV.value;
     } 
   }
      
   // now, we can take smallest positive distance.
   
   // Initialize
   //
   G4double      distance = kInfinity;   

   // Find intersections and choose nearest one
   //
   G4VTwistSurface *surfaces[6];

   surfaces[0] = fSide0;
   surfaces[1] = fSide90 ;
   surfaces[2] = fSide180 ;
   surfaces[3] = fSide270 ;
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;
   
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   G4int i;
   for (i=0; i < 6 ; i++)
   {
#ifdef G4TWISTDEBUG
      G4cout << G4endl << "surface " << i << ": " << G4endl << G4endl ;
#endif
      G4double tmpdistance = surfaces[i]->DistanceToIn(p, v, xx);
#ifdef G4TWISTDEBUG
      G4cout << "Solid DistanceToIn : distance = " << tmpdistance << G4endl ; 
      G4cout << "intersection point = " << xx << G4endl ;
#endif 
      if (tmpdistance < distance)
      {
         distance = tmpdistance;
         bestxx = xx;
      }
   }

#ifdef G4TWISTDEBUG
   G4cout << "best distance = " << distance << G4endl ;
#endif

   *tmpdist = distance;
   // timer.Stop();
   return fLastDistanceToInWithV.value;
}


//=====================================================================
//* DistanceToIn (p) --------------------------------------------------

G4double G4VTwistedFaceted::DistanceToIn (const G4ThreeVector& p) const
{
   // DistanceToIn(p):
   // Calculate distance to surface of shape from `outside',
   // allowing for tolerance
   //
   
   //
   // checking last value
   //
   
   G4ThreeVector *tmpp;
   G4double      *tmpdist;
   if (fLastDistanceToIn.p == p)
   {
      return fLastDistanceToIn.value;
   }
   else
   {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToIn.p));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToIn.value));
      tmpp->set(p.x(), p.y(), p.z());
   }

   //
   // Calculate DistanceToIn(p) 
   //
   
   EInside currentside = Inside(p);

   switch (currentside)
   {
      case (kInside) :
      {
      }

      case (kSurface) :
      {
         *tmpdist = 0.;
         return fLastDistanceToIn.value;
      }

      case (kOutside) :
      {
         // Initialize
         //
         G4double      distance = kInfinity;   

         // Find intersections and choose nearest one
         //
         G4VTwistSurface *surfaces[6];

         surfaces[0] = fSide0;
         surfaces[1] = fSide90 ;
         surfaces[2] = fSide180 ;
         surfaces[3] = fSide270 ;
         surfaces[4] = fLowerEndcap;
         surfaces[5] = fUpperEndcap;

         G4int i;
         G4ThreeVector xx;
         G4ThreeVector bestxx;
         for (i=0; i< 6; i++)
         {
            G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
            if (tmpdistance < distance)
            {
               distance = tmpdistance;
               bestxx = xx;
            }
         }
         *tmpdist = distance;
         return fLastDistanceToIn.value;
      }

      default :
      {
         G4Exception("G4VTwistedFaceted::DistanceToIn(p)", "GeomSolids0003",
                     FatalException, "Unknown point location!");
      }
   } // switch end

   return 0;
}


//=====================================================================
//* DistanceToOut (p, v) ----------------------------------------------

G4double
G4VTwistedFaceted::DistanceToOut( const G4ThreeVector& p,
                                  const G4ThreeVector& v,
                                  const G4bool calcNorm,
                                        G4bool *validNorm,
                                        G4ThreeVector *norm ) const
{
   // DistanceToOut (p, v):
   // Calculate distance to surface of shape from `inside'
   // along with the v, allowing for tolerance.
   // The function returns kInfinity if no intersection or
   // just grazing within tolerance.

   //
   // checking last value
   //
   
   G4ThreeVector *tmpp;
   G4ThreeVector *tmpv;
   G4double      *tmpdist;
   if (fLastDistanceToOutWithV.p == p && fLastDistanceToOutWithV.vec == v  )
   {
      return fLastDistanceToOutWithV.value;
   }
   else
   {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToOutWithV.p));
      tmpv    = const_cast<G4ThreeVector*>(&(fLastDistanceToOutWithV.vec));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToOutWithV.value));
      tmpp->set(p.x(), p.y(), p.z());
      tmpv->set(v.x(), v.y(), v.z());
   }

   //
   // Calculate DistanceToOut(p,v)
   //
   
   EInside currentside = Inside(p);

   if (currentside == kOutside)
   {
   }
   else if (currentside == kSurface)
   {
      // particle is just on a boundary.
      // if the particle is exiting from the volume, return 0
      //
      G4ThreeVector normal = SurfaceNormal(p);
      G4VTwistSurface *blockedsurface = fLastNormal.surface[0];
      if (normal*v > 0)
      {
            if (calcNorm)
            {
               *norm = (blockedsurface->GetNormal(p, true));
               *validNorm = blockedsurface->IsValidNorm();
            }
            *tmpdist = 0.;
            // timer.Stop();
            return fLastDistanceToOutWithV.value;
      }
   }
      
   // now, we can take smallest positive distance.
   
   // Initialize
   G4double      distance = kInfinity;
       
   // find intersections and choose nearest one.
   G4VTwistSurface *surfaces[6];

   surfaces[0] = fSide0;
   surfaces[1] = fSide90 ;
   surfaces[2] = fSide180 ;
   surfaces[3] = fSide270 ;
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;

   G4int i;
   G4int besti = -1;
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   for (i=0; i< 6 ; i++) {
      G4double tmpdistance = surfaces[i]->DistanceToOut(p, v, xx);
      if (tmpdistance < distance)
      {
         distance = tmpdistance;
         bestxx = xx; 
         besti = i;
      }
   }

   if (calcNorm)
   {
      if (besti != -1)
      {
         *norm = (surfaces[besti]->GetNormal(p, true));
         *validNorm = surfaces[besti]->IsValidNorm();
      }
   }

   *tmpdist = distance;
   // timer.Stop();
   return fLastDistanceToOutWithV.value;
}


//=====================================================================
//* DistanceToOut (p) -------------------------------------------------

G4double G4VTwistedFaceted::DistanceToOut( const G4ThreeVector& p ) const
{
   // DistanceToOut(p):
   // Calculate distance to surface of shape from `inside', 
   // allowing for tolerance

   //
   // checking last value
   //

   G4ThreeVector *tmpp;
   G4double      *tmpdist;

   if (fLastDistanceToOut.p == p)
   {
      return fLastDistanceToOut.value;
   }
   else
   {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToOut.p));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToOut.value));
      tmpp->set(p.x(), p.y(), p.z());
   }
   
   //
   // Calculate DistanceToOut(p)
   //
   
   EInside currentside = Inside(p);
   G4double     retval = kInfinity;   

   switch (currentside)
   {
      case (kOutside) :
      {
#ifdef G4SPECSDEBUG
        G4int oldprc = G4cout.precision(16) ;
        G4cout << G4endl ;
        DumpInfo();
        G4cout << "Position:"  << G4endl << G4endl ;
        G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl ;
        G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl ;
        G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl ;
        G4cout.precision(oldprc) ;
        G4Exception("G4VTwistedFaceted::DistanceToOut(p)", "GeomSolids1002",
                    JustWarning, "Point p is outside !?" );
#endif
        break;
      }
      case (kSurface) :
      {
        *tmpdist = 0.;
        retval = fLastDistanceToOut.value;
        break;
      }
      
      case (kInside) :
      {
        // Initialize
        //
        G4double      distance = kInfinity;
   
        // find intersections and choose nearest one
        //
        G4VTwistSurface *surfaces[6];

        surfaces[0] = fSide0;
        surfaces[1] = fSide90 ;
        surfaces[2] = fSide180 ;
        surfaces[3] = fSide270 ;
        surfaces[4] = fLowerEndcap;
        surfaces[5] = fUpperEndcap;

        G4int i;
        G4ThreeVector xx;
        G4ThreeVector bestxx;
        for (i=0; i< 6; i++)
        {
          G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
          if (tmpdistance < distance)
          {
            distance = tmpdistance;
            bestxx = xx;
          }
        }
        *tmpdist = distance;
   
        retval = fLastDistanceToOut.value;
        break;
      }
      
      default :
      {
        G4Exception("G4VTwistedFaceted::DistanceToOut(p)", "GeomSolids0003",
                    FatalException, "Unknown point location!");
        break;
      }
   } // switch end

   return retval;
}


//=====================================================================
//* StreamInfo --------------------------------------------------------

std::ostream& G4VTwistedFaceted::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4VTwistedFaceted\n"
     << " Parameters: \n"
     << "  polar angle theta = "   <<  fTheta/degree      << " deg" << G4endl
     << "  azimuthal angle phi = "  << fPhi/degree        << " deg" << G4endl  
     << "  tilt angle  alpha = "   << fAlph/degree        << " deg" << G4endl  
     << "  TWIST angle = "         << fPhiTwist/degree    << " deg" << G4endl  
     << "  Half length along y (lower endcap) = "         << fDy1/cm << " cm"
     << G4endl 
     << "  Half length along x (lower endcap, bottom) = " << fDx1/cm << " cm"
     << G4endl 
     << "  Half length along x (lower endcap, top) = "    << fDx2/cm << " cm"
     << G4endl 
     << "  Half length along y (upper endcap) = "         << fDy2/cm << " cm"
     << G4endl 
     << "  Half length along x (upper endcap, bottom) = " << fDx3/cm << " cm"
     << G4endl 
     << "  Half length along x (upper endcap, top) = "    << fDx4/cm << " cm"
     << G4endl 
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


//=====================================================================
//* DiscribeYourselfTo ------------------------------------------------

void G4VTwistedFaceted::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddSolid (*this);
}


//=====================================================================
//* GetExtent ---------------------------------------------------------

G4VisExtent G4VTwistedFaceted::GetExtent() const 
{
  G4double maxRad = std::sqrt( fDx*fDx + fDy*fDy);

  return G4VisExtent(-maxRad, maxRad ,
                     -maxRad, maxRad ,
                     -fDz, fDz );
}


//=====================================================================
//* CreateSurfaces ----------------------------------------------------

void G4VTwistedFaceted::CreateSurfaces() 
{
   
  // create 6 surfaces of TwistedTub.

  if ( fDx1 == fDx2 && fDx3 == fDx4 )    // special case : Box
  {
    fSide0   = new G4TwistBoxSide("0deg",   fPhiTwist, fDz, fTheta, fPhi,
                                        fDy1, fDx1, fDx1, fDy2, fDx3, fDx3, fAlph, 0.*deg);
    fSide180 = new G4TwistBoxSide("180deg", fPhiTwist, fDz, fTheta, fPhi+pi,
                                        fDy1, fDx1, fDx1, fDy2, fDx3, fDx3, fAlph, 180.*deg);
  }
  else   // default general case
  {
    fSide0   = new G4TwistTrapAlphaSide("0deg"   ,fPhiTwist, fDz, fTheta,
                      fPhi, fDy1, fDx1, fDx2, fDy2, fDx3, fDx4, fAlph, 0.*deg);
    fSide180 = new G4TwistTrapAlphaSide("180deg", fPhiTwist, fDz, fTheta,
                 fPhi+pi, fDy1, fDx2, fDx1, fDy2, fDx4, fDx3, fAlph, 180.*deg);
  }

  // create parallel sides
  //
  fSide90 = new G4TwistTrapParallelSide("90deg",  fPhiTwist, fDz, fTheta,
                      fPhi, fDy1, fDx1, fDx2, fDy2, fDx3, fDx4, fAlph, 0.*deg);
  fSide270 = new G4TwistTrapParallelSide("270deg", fPhiTwist, fDz, fTheta,
                 fPhi+pi, fDy1, fDx2, fDx1, fDy2, fDx4, fDx3, fAlph, 180.*deg);

  // create endcaps
  //
  fUpperEndcap = new G4TwistTrapFlatSide("UpperCap",fPhiTwist, fDx3, fDx4, fDy2,
                                    fDz, fAlph, fPhi, fTheta,  1 );
  fLowerEndcap = new G4TwistTrapFlatSide("LowerCap",fPhiTwist, fDx1, fDx2, fDy1,
                                    fDz, fAlph, fPhi, fTheta, -1 );
 
  // Set neighbour surfaces
  
  fSide0->SetNeighbours(  fSide270 , fLowerEndcap , fSide90  , fUpperEndcap );
  fSide90->SetNeighbours( fSide0   , fLowerEndcap , fSide180 , fUpperEndcap );
  fSide180->SetNeighbours(fSide90  , fLowerEndcap , fSide270 , fUpperEndcap );
  fSide270->SetNeighbours(fSide180 , fLowerEndcap , fSide0   , fUpperEndcap );
  fUpperEndcap->SetNeighbours( fSide180, fSide270 , fSide0 , fSide90  );
  fLowerEndcap->SetNeighbours( fSide180, fSide270 , fSide0 , fSide90  );

}


//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4VTwistedFaceted::GetEntityType() const
{
  return G4String("G4VTwistedFaceted");
}


//=====================================================================
//* GetPolyhedron -----------------------------------------------------

G4Polyhedron* G4VTwistedFaceted::GetPolyhedron() const
{
  if (!fpPolyhedron ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      G4AutoLock l(&polyhedronMutex);
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
      fRebuildPolyhedron = false;
      l.unlock();
    }

  return fpPolyhedron;
}


//=====================================================================
//* GetPointInSolid ---------------------------------------------------

G4ThreeVector G4VTwistedFaceted::GetPointInSolid(G4double z) const
{


  // this routine is only used for a test
  // can be deleted ...

  if ( z == fDz ) z -= 0.1*fDz ;
  if ( z == -fDz ) z += 0.1*fDz ;

  G4double phi = z/(2*fDz)*fPhiTwist ;

  return G4ThreeVector(fdeltaX * phi/fPhiTwist, fdeltaY * phi/fPhiTwist, z ) ;
}


//=====================================================================
//* GetPointOnSurface -------------------------------------------------

G4ThreeVector G4VTwistedFaceted::GetPointOnSurface() const
{

  G4double  phi = G4RandFlat::shoot(-fPhiTwist/2.,fPhiTwist/2.);
  G4double u , umin, umax ;  //  variable for twisted surfaces
  G4double y  ;              //  variable for flat surface (top and bottom)

  // Compute the areas. Attention: Only correct for trapezoids
  // where the twisting is done along the z-axis. In the general case
  // the computed surface area is more difficult. However this simplification
  // does not affect the tracking through the solid. 
 
  G4double a1   = fSide0->GetSurfaceArea();
  G4double a2   = fSide90->GetSurfaceArea();
  G4double a3   = fSide180->GetSurfaceArea() ;
  G4double a4   = fSide270->GetSurfaceArea() ;
  G4double a5   = fLowerEndcap->GetSurfaceArea() ;
  G4double a6   = fUpperEndcap->GetSurfaceArea() ;

#ifdef G4TWISTDEBUG
  G4cout << "Surface 0   deg = " << a1 << G4endl ;
  G4cout << "Surface 90  deg = " << a2 << G4endl ;
  G4cout << "Surface 180 deg = " << a3 << G4endl ;
  G4cout << "Surface 270 deg = " << a4 << G4endl ;
  G4cout << "Surface Lower   = " << a5 << G4endl ;
  G4cout << "Surface Upper   = " << a6 << G4endl ;
#endif 

  G4double chose = G4RandFlat::shoot(0.,a1 + a2 + a3 + a4 + a5 + a6) ;

  if(chose < a1)
  {

    umin = fSide0->GetBoundaryMin(phi) ;
    umax = fSide0->GetBoundaryMax(phi) ;
    u = G4RandFlat::shoot(umin,umax) ;

    return  fSide0->SurfacePoint(phi, u, true) ;   // point on 0deg surface
  }

  else if( (chose >= a1) && (chose < a1 + a2 ) )
  {

    umin = fSide90->GetBoundaryMin(phi) ;
    umax = fSide90->GetBoundaryMax(phi) ;
    
    u = G4RandFlat::shoot(umin,umax) ;

    return fSide90->SurfacePoint(phi, u, true);   // point on 90deg surface
  }

  else if( (chose >= a1 + a2 ) && (chose < a1 + a2 + a3 ) )
  {

    umin = fSide180->GetBoundaryMin(phi) ;
    umax = fSide180->GetBoundaryMax(phi) ;
    u = G4RandFlat::shoot(umin,umax) ;

     return fSide180->SurfacePoint(phi, u, true); // point on 180 deg surface
  }

  else if( (chose >= a1 + a2 + a3  ) && (chose < a1 + a2 + a3 + a4  ) )
  {

    umin = fSide270->GetBoundaryMin(phi) ;
    umax = fSide270->GetBoundaryMax(phi) ;
    u = G4RandFlat::shoot(umin,umax) ;

    return fSide270->SurfacePoint(phi, u, true); // point on 270 deg surface
  }

  else if( (chose >= a1 + a2 + a3 + a4  ) && (chose < a1 + a2 + a3 + a4 + a5 ) )
  {

    y = G4RandFlat::shoot(-fDy1,fDy1) ;
    umin = fLowerEndcap->GetBoundaryMin(y) ;
    umax = fLowerEndcap->GetBoundaryMax(y) ;
    u = G4RandFlat::shoot(umin,umax) ;

    return fLowerEndcap->SurfacePoint(u,y,true); // point on lower endcap
  }
  else {

    y = G4RandFlat::shoot(-fDy2,fDy2) ;
    umin = fUpperEndcap->GetBoundaryMin(y) ;
    umax = fUpperEndcap->GetBoundaryMax(y) ;
    u = G4RandFlat::shoot(umin,umax) ;

    return fUpperEndcap->SurfacePoint(u,y,true) ; // point on upper endcap

  }
}


//=====================================================================
//* CreatePolyhedron --------------------------------------------------

G4Polyhedron* G4VTwistedFaceted::CreatePolyhedron () const 
{
  // number of meshes
  const G4int k =
    G4int(G4Polyhedron::GetNumberOfRotationSteps() * fPhiTwist / twopi) + 2;
  const G4int n = k;

  const G4int nnodes = 4*(k-1)*(n-2) + 2*k*k ;
  const G4int nfaces = 4*(k-1)*(n-1) + 2*(k-1)*(k-1) ;

  G4Polyhedron *ph=new G4Polyhedron;
  typedef G4double G4double3[3];
  typedef G4int G4int4[4];
  G4double3* xyz = new G4double3[nnodes];  // number of nodes 
  G4int4*  faces = new G4int4[nfaces] ;    // number of faces

  fLowerEndcap->GetFacets(k,k,xyz,faces,0) ;
  fUpperEndcap->GetFacets(k,k,xyz,faces,1) ;
  fSide270->GetFacets(k,n,xyz,faces,2) ;
  fSide0->GetFacets(k,n,xyz,faces,3) ;
  fSide90->GetFacets(k,n,xyz,faces,4) ;
  fSide180->GetFacets(k,n,xyz,faces,5) ;

  ph->createPolyhedron(nnodes,nfaces,xyz,faces);

  return ph;
}
