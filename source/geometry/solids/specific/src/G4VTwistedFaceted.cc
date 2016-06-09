//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4VTwistedFaceted.cc,v 1.2 2005/04/04 11:56:59 gcosmo Exp $
// GEANT4 tag $Name: geant4-07-01 $
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

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4SolidExtentList.hh"
#include "G4ClippablePolygon.hh"
#include "G4VPVParameterisation.hh"
#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"
#include "G4NURBS.hh"
#include "G4NURBStube.hh"
#include "G4NURBScylinder.hh"
#include "G4NURBStubesector.hh"

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
  : G4VSolid(pname), 
    fLowerEndcap(0), fUpperEndcap(0), fSide0(0),
    fSide90(0), fSide180(0), fSide270(0),
    fpPolyhedron(0)
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
      G4cerr << "ERROR - G4VTwistedFaceted::G4VTwistedFaceted(): "
             << GetName() << G4endl
             << "        Not planar ! - " << G4endl 
             << "fDy2 is " << fDy2 << " but should be "
             << pDytmp << "." << G4endl ;
      G4Exception("G4VTwistedFaceted::G4VTwistedFaceted()", "InvalidSetup",
                  FatalException, "Not planar surface in untwisted Trapezoid.");
    }
  }

#ifdef G4SPECSDEBUG
  if ( fDx1 == fDx2 && fDx3 == fDx4 )
  { 
      G4cout << "Trapezoid is a box" << G4endl ;
  }
  
#endif

  if ( (  fDx1 == fDx2 && fDx3 != fDx4 ) || ( fDx1 != fDx2 && fDx3 == fDx4 ) )
  {
    G4cerr << "ERROR - G4VTwistedFaceted::G4VTwistedFaceted(): "
           << GetName() << G4endl
           << "        Not planar ! - " << G4endl 
           << "One endcap is rectengular, the other is a trapezoid." << G4endl
           << "For planarity reasons they have to be rectangles or trapezoids "
           << G4endl
           << "on both sides."
           << G4endl ;
    G4Exception("G4VTwistedFaceted::G4VTwistedFaceted()", "InvalidSetup",
                FatalException, "Not planar surface in untwisted Trapezoid.");
  }

  // twist angle
  //
  fPhiTwist = PhiTwist ;

  // tilt angle
  //
  fAlph = - pAlph ;  // important: definition of angle alpha
                     //            is different in equations !!!
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
         && ( fDy1   > 2*kCarTolerance)
         && ( fDy2   > 2*kCarTolerance)
         && ( fDz   > 2*kCarTolerance) 
         && ( std::fabs(fPhiTwist) > 2*kAngTolerance )
         && ( std::fabs(fPhiTwist) < pi/2 )
         && ( std::fabs(fAlph) < pi/2 )
         && ( fTheta < pi/2 && fTheta >= 0 ) )
      )
  {
    G4cerr << "ERROR - G4VTwistedFaceted()::G4VTwistedFaceted(): "
           << GetName() << G4endl
           << "        Dimensions too small or too big! - " << G4endl 
           << "fDx 1-4 = " << fDx1/cm << ", " << fDx2/cm << ", "
           << fDx3/cm << ", " << fDx4/cm << " cm" << G4endl 
           << "fDy 1-2 = " << fDy1/cm << ", " << fDy2/cm << ", "
           << " cm" << G4endl 
           << "fDz = " << fDz/cm << " cm" << G4endl 
           << " twistangle " << fPhiTwist/deg << " deg" << G4endl 
           << " phi,theta = " << fPhi/deg << ", "  << fTheta/deg
           << " deg" << G4endl ;
    G4Exception("G4TwistedTrap::G4VTwistedFaceted()",
                "InvalidSetup", FatalException,
                "Invalid dimensions. Too small, or twist angle too big.");
  }
  CreateSurfaces();
  fCubicVolume = 2 * fDz * ( ( fDx1 + fDx2 ) * fDy1 + ( fDx3 + fDx4 ) * fDy2 );
}


//=====================================================================
//* destructor --------------------------------------------------------

G4VTwistedFaceted::~G4VTwistedFaceted()
{
  if (fLowerEndcap) delete fLowerEndcap ;
  if (fUpperEndcap) delete fUpperEndcap ;

  if (fSide0)       delete fSide0 ;
  if (fSide90)      delete fSide90 ;
  if (fSide180)     delete fSide180 ;
  if (fSide270)     delete fSide270 ;
  if (fpPolyhedron) delete fpPolyhedron;
}

//=====================================================================
//* ComputeDimensions -------------------------------------------------

void G4VTwistedFaceted::ComputeDimensions(G4VPVParameterisation* ,
                                          const G4int ,
                                          const G4VPhysicalVolume* )
{
  G4Exception("G4VTwistedFaceted::ComputeDimensions()",
              "NotSupported", FatalException,
              "G4VTwistedFaceted does not support Parameterisation.");
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
  G4double maxRad = std::sqrt( fDx*fDx + fDy*fDy);

  if (!pTransform.IsRotated())
    {
      // Special case handling for unrotated boxes
      // Compute x/y/z mins and maxs respecting limits, with early returns
      // if outside limits. Then switch() on pAxis
      
      G4double xoffset,xMin,xMax;
      G4double yoffset,yMin,yMax;
      G4double zoffset,zMin,zMax;

      xoffset = pTransform.NetTranslation().x() ;
      xMin    = xoffset - maxRad ;
      xMax    = xoffset + maxRad ;

      if (pVoxelLimit.IsXLimited())
        {
          if ( xMin > pVoxelLimit.GetMaxXExtent()+kCarTolerance || 
               xMax < pVoxelLimit.GetMinXExtent()-kCarTolerance  ) return false;
          else
            {
              if (xMin < pVoxelLimit.GetMinXExtent())
                {
                  xMin = pVoxelLimit.GetMinXExtent() ;
                }
              if (xMax > pVoxelLimit.GetMaxXExtent())
                {
                  xMax = pVoxelLimit.GetMaxXExtent() ;
                }
            }
        }
      yoffset = pTransform.NetTranslation().y() ;
      yMin    = yoffset - maxRad ;
      yMax    = yoffset + maxRad ;
      
      if (pVoxelLimit.IsYLimited())
        {
          if ( yMin > pVoxelLimit.GetMaxYExtent()+kCarTolerance ||
               yMax < pVoxelLimit.GetMinYExtent()-kCarTolerance  ) return false;
          else
            {
              if (yMin < pVoxelLimit.GetMinYExtent())
                {
                  yMin = pVoxelLimit.GetMinYExtent() ;
                }
              if (yMax > pVoxelLimit.GetMaxYExtent())
                {
                  yMax = pVoxelLimit.GetMaxYExtent() ;
                }
            }
        }
      zoffset = pTransform.NetTranslation().z() ;
      zMin    = zoffset - fDz ;
      zMax    = zoffset + fDz ;
      
      if (pVoxelLimit.IsZLimited())
        {
          if ( zMin > pVoxelLimit.GetMaxZExtent()+kCarTolerance ||
               zMax < pVoxelLimit.GetMinZExtent()-kCarTolerance  ) return false;
          else
            {
              if (zMin < pVoxelLimit.GetMinZExtent())
                {
                  zMin = pVoxelLimit.GetMinZExtent() ;
                }
              if (zMax > pVoxelLimit.GetMaxZExtent())
                {
                  zMax = pVoxelLimit.GetMaxZExtent() ;
                }
            }
        }
      switch (pAxis)
        {
        case kXAxis:
          pMin = xMin ;
          pMax = xMax ;
          break ;
        case kYAxis:
          pMin=yMin;
          pMax=yMax;
          break;
        case kZAxis:
        pMin=zMin;
        pMax=zMax;
        break;
        default:
          break;
        }
      pMin -= kCarTolerance ;
      pMax += kCarTolerance ;
      
      return true;
  }
  else  // General rotated case - create and clip mesh to boundaries
    {
      G4bool existsAfterClip = false ;
      G4ThreeVectorList* vertices ;

      pMin = +kInfinity ;
      pMax = -kInfinity ;
    
    // Calculate rotated vertex coordinates
      
      vertices = CreateRotatedVertices(pTransform) ;
      ClipCrossSection(vertices,0,pVoxelLimit,pAxis,pMin,pMax) ;
      ClipCrossSection(vertices,4,pVoxelLimit,pAxis,pMin,pMax) ;
      ClipBetweenSections(vertices,0,pVoxelLimit,pAxis,pMin,pMax) ;
      
      if (pVoxelLimit.IsLimited(pAxis) == false) 
        {  
          if ( pMin != kInfinity || pMax != -kInfinity ) 
            {
              existsAfterClip = true ;

              // Add 2*tolerance to avoid precision troubles

              pMin           -= kCarTolerance;
              pMax           += kCarTolerance;
            }
        }      
      else
        {
          G4ThreeVector clipCentre(
            ( pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
            ( pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
            ( pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5 );
          
      if ( pMin != kInfinity || pMax != -kInfinity )
        {
          existsAfterClip = true ;
          
          // Check to see if endpoints are in the solid
          
          clipCentre(pAxis) = pVoxelLimit.GetMinExtent(pAxis);
          
          if (Inside(pTransform.Inverse().TransformPoint(clipCentre))
              != kOutside)
            {
              pMin = pVoxelLimit.GetMinExtent(pAxis);
            }
          else
            {
              pMin -= kCarTolerance;
            }
          clipCentre(pAxis) = pVoxelLimit.GetMaxExtent(pAxis);
          
          if (Inside(pTransform.Inverse().TransformPoint(clipCentre))
              != kOutside)
            {
              pMax = pVoxelLimit.GetMaxExtent(pAxis);
            }
          else
            {
              pMax += kCarTolerance;
            }
        }
      
      // Check for case where completely enveloping clipping volume
      // If point inside then we are confident that the solid completely
      // envelopes the clipping volume. Hence set min/max extents according
      // to clipping volume extents along the specified axis.
      
      else if (Inside(pTransform.Inverse().TransformPoint(clipCentre))
               != kOutside)
        {
          existsAfterClip = true ;
          pMin            = pVoxelLimit.GetMinExtent(pAxis) ;
          pMax            = pVoxelLimit.GetMaxExtent(pAxis) ;
        }
        } 
      delete vertices;
      return existsAfterClip;
    } 
  
  
}

G4ThreeVectorList* G4VTwistedFaceted::
CreateRotatedVertices(const G4AffineTransform& pTransform) const
{

  G4ThreeVectorList* vertices = new G4ThreeVectorList();
  vertices->reserve(8);

  if (vertices)
  {

    G4double maxRad = std::sqrt( fDx*fDx + fDy*fDy);

    G4ThreeVector vertex0(-maxRad,-maxRad,-fDz) ;
    G4ThreeVector vertex1(maxRad,-maxRad,-fDz) ;
    G4ThreeVector vertex2(maxRad,maxRad,-fDz) ;
    G4ThreeVector vertex3(-maxRad,maxRad,-fDz) ;
    G4ThreeVector vertex4(-maxRad,-maxRad,fDz) ;
    G4ThreeVector vertex5(maxRad,-maxRad,fDz) ;
    G4ThreeVector vertex6(maxRad,maxRad,fDz) ;
    G4ThreeVector vertex7(-maxRad,maxRad,fDz) ;

    vertices->push_back(pTransform.TransformPoint(vertex0));
    vertices->push_back(pTransform.TransformPoint(vertex1));
    vertices->push_back(pTransform.TransformPoint(vertex2));
    vertices->push_back(pTransform.TransformPoint(vertex3));
    vertices->push_back(pTransform.TransformPoint(vertex4));
    vertices->push_back(pTransform.TransformPoint(vertex5));
    vertices->push_back(pTransform.TransformPoint(vertex6));
    vertices->push_back(pTransform.TransformPoint(vertex7));
  }
  else
  {
    DumpInfo();
    G4Exception("G4VTwistedFaceted::CreateRotatedVertices()",
                "FatalError", FatalException,
                "Error in allocation of vertices. Out of memory !");
  }
  return vertices;
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
      tmpin = const_cast<EInside*>(&(fLastInside.inside));
      tmpp->set(p.x(), p.y(), p.z());
   }

   *tmpin = kOutside ;

   G4double phi = p.z()/(2*fDz) * fPhiTwist ;  // rotate the point to z=0
   G4double cphi = std::cos(-phi) ;
   G4double sphi = std::sin(-phi) ;
   G4double posx = p.x() * cphi - p.y() * sphi   ;
   G4double posy = p.x() * sphi + p.y() * cphi   ;
   G4double posz = p.z()  ;

   posx += fdeltaX * ( -phi/fPhiTwist) ;  // shift
   posy += fdeltaY * ( -phi/fPhiTwist) ;

   G4double xMin = Xcoef(posy,phi,fTAlph) - 2*Xcoef(posy,phi,0.) ; 
   G4double xMax = Xcoef(posy,phi,fTAlph) ;  

   G4double yMax = GetValueB(phi)/2. ;  // b(phi)/2 is limit
   G4double yMin = -yMax ;

#ifdef G4SPECSDEBUG

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

#ifdef G4SPECSDEBUG
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
   G4VSurface    **tmpsurface = const_cast<G4VSurface**>(fLastNormal.surface);
   tmpp->set(p.x(), p.y(), p.z());

   G4double      distance = kInfinity;

   G4VSurface *surfaces[6];

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
   G4VSurface *surfaces[6];

   surfaces[0] = fSide0;
   surfaces[1] = fSide90 ;
   surfaces[2] = fSide180 ;
   surfaces[3] = fSide270 ;
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;
   
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   G4int i;
   G4int besti = -1;
   for (i=0; i < 6 ; i++)
   {

#ifdef G4SPECSDEBUG
      G4cout << G4endl << "surface " << i << ": " << G4endl << G4endl ;
#endif
      G4double tmpdistance = surfaces[i]->DistanceToIn(p, v, xx);
#ifdef G4SPECSDEBUG
      G4cout << "Solid DistanceToIn : distance = " << tmpdistance << G4endl ; 
      G4cout << "intersection point = " << xx << G4endl ;
#endif 
      if (tmpdistance < distance)
      {
         distance = tmpdistance;
         bestxx = xx;
         besti = i;
      }
   }

#ifdef G4SPECSDEBUG
   G4cout << "best distance = " << distance << G4endl ;
#endif

   *tmpdist = distance;
   // timer.Stop();
   return fLastDistanceToInWithV.value;
}


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
         G4VSurface *surfaces[6];

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
         *tmpdist = distance;
         return fLastDistanceToIn.value;
      }

      default :
      {
         G4Exception("G4VTwistedFaceted::DistanceToIn(p)", "InvalidCondition",
                     FatalException, "Unknown point location!");
      }
   } // switch end

   return kInfinity;
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
      G4VSurface *blockedsurface = fLastNormal.surface[0];
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
   G4VSurface *surfaces[6];

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

   switch (currentside)
   {
      case (kOutside) :
      {
      }
      case (kSurface) :
      {
        *tmpdist = 0.;
         return fLastDistanceToOut.value;
      }
      
      case (kInside) :
      {
         // Initialize
         //
         G4double      distance = kInfinity;
   
         // find intersections and choose nearest one
         //
         G4VSurface *surfaces[6];

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
         *tmpdist = distance;
   
         return fLastDistanceToOut.value;
      }
      
      default :
      {
         G4Exception("G4VTwistedFaceted::DistanceToOut(p)", "InvalidCondition",
                     FatalException, "Unknown point location!");
      }
   } // switch end

   return 0;
}


//=====================================================================
//* StreamInfo --------------------------------------------------------

std::ostream& G4VTwistedFaceted::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
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
//* CreateNUBS --------------------------------------------------------

G4NURBS* G4VTwistedFaceted::CreateNURBS () const 
{
  G4double maxRad = std::sqrt( fDx*fDx + fDy*fDy);

  return new G4NURBStube(maxRad, maxRad, fDz); 
   // Tube for now!!!
}


//=====================================================================
//* CreateSurfaces ----------------------------------------------------

void G4VTwistedFaceted::CreateSurfaces() 
{
   
  // create 6 surfaces of TwistedTub.

  if ( fDx1 == fDx2 && fDx3 == fDx4 )    // special case : Box
  {
    fSide0   = new G4TwistedTrapBoxSide("0deg",   fPhiTwist, fDz, fTheta, fPhi,
                                        fDy1, fDx1, fDy2, fDx3, fAlph, 0.*deg);
    fSide180 = new G4TwistedTrapBoxSide("180deg", fPhiTwist, fDz, fTheta,
                            fPhi+pi , fDy1, fDx1, fDy2, fDx3, fAlph, 180.*deg);
  }
  else   // default general case
  {
    fSide0   = new G4TwistedTrapAlphaSide("0deg"   ,fPhiTwist, fDz, fTheta,
                      fPhi, fDy1, fDx1, fDx2, fDy2, fDx3, fDx4, fAlph, 0.*deg);
    fSide180 = new G4TwistedTrapAlphaSide("180deg", fPhiTwist, fDz, fTheta,
                 fPhi+pi, fDy1, fDx2, fDx1, fDy2, fDx4, fDx3, fAlph, 180.*deg);
  }

  // create parallel sides
  //
  fSide90 = new G4TwistedTrapParallelSide("90deg",  fPhiTwist, fDz, fTheta,
                      fPhi, fDy1, fDx1, fDx2, fDy2, fDx3, fDx4, fAlph, 0.*deg);
  fSide270 = new G4TwistedTrapParallelSide("270deg", fPhiTwist, fDz, fTheta,
                 fPhi+pi, fDy1, fDx2, fDx1, fDy2, fDx4, fDx3, fAlph, 180.*deg);

  // create endcaps
  //
  fUpperEndcap = new G4FlatTrapSide("UpperCap",fPhiTwist, fDx3, fDx4, fDy2,
                                    fDz, fAlph, fPhi, fTheta,  1 );
  fLowerEndcap = new G4FlatTrapSide("LowerCap",fPhiTwist, fDx1, fDx2, fDy1,
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
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      delete fpPolyhedron;
      fpPolyhedron = CreatePolyhedron();
    }
  return fpPolyhedron;
}
