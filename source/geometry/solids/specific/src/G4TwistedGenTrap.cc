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
// $Id: G4TwistedGenTrap.cc,v 1.1 2005-03-02 15:40:39 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedGenTrap.cc
//
// Author:
//
//   04-Nov-2004 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

// #define DISTANCETOIN

#include "G4TwistedGenTrap.hh"

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

G4TwistedGenTrap::G4TwistedGenTrap(const G4String &pname,         // Name of instance
				   G4double      PhiTwist,    // twist angle
				   G4double      pDz,         // half z lenght
				   G4double      pTheta,      // direction between end planes
				   G4double      pPhi,        // defined by polar and azimutal angles.
				   G4double      pDy1,        // half y length at -pDz
				   G4double      pDx1,        // half x length at -pDz,-pDy
				   G4double      pDx2,        // half x length at -pDz,+pDy
				   G4double      pAlph1,      // tilt angle at -pDz
				   G4double      pDy2,        // half y length at +pDz
				   G4double      pDx3,        // half x length at +pDz,-pDy
				   G4double      pDx4,        // half x length at +pDz,+pDy
				   G4double      pAlph2      // tilt angle at +pDzconst G4String &pname,
				   )
  : G4VSolid(pname), 
     fLowerEndcap(0), fUpperEndcap(0), fSide0(0),
    fSide90(0), fSide180(0), fSide270(0)
{

  G4double pDytmp ;
  fDx1 = pDx1 ;
  fDx2 = pDx2 ;
  fDx3 = pDx3 ;
  fDx4 = pDx4 ;
  fDz   = pDz ;
  
  // planarity check ( or in the constructor of the solid ?? )
  if ( fDx1 != fDx2 && fDx3 != fDx4 ) {
    pDytmp = fDy1 * ( fDx3 - fDx4 ) / ( fDx1 - fDx2 ) ;
    if ( std::fabs(pDytmp - fDy2) > kCarTolerance ) {
      G4cerr << "ERROR - G4TwistedGenTrap::G4TwistedGenTrap(): " << GetName() << G4endl
             << "        Not planar ! - " << G4endl 
             << "fDy2 is " << fDy2 << " but should be " << pDytmp << "." << G4endl ;
      
      G4Exception("G4TwistedGenTrap::G4TwistedGenTrap()", "InvalidSetup",
                  FatalException, "Not planar surface in untwisted Trapezoid.");
    }
    else {
      G4cerr << "ERROR - G4TwistedGenTrap::G4TwistedGenTrap(): " << GetName() << G4endl
             << "        Not implemented ! - " << G4endl 
             << "Use G4TwistedBox instead." << G4endl ;
      
      G4Exception("G4TwistedGenTrap::G4TwistedGenTrap()", "InvalidSetup",
                  FatalException, "Not implemented ! - Use G4TwistedBox instead.");
      
    }
  }
    
  fPhiTwist = PhiTwist ;
  
  fAlph = pAlph1 ;
  // add test pAlph1 == pAlph2 !!!!
  fTAlph = std::tan(fAlph) ;
  
  fTheta = pTheta ;
  fPhi   = pPhi ;


  CreateSurfaces();
  fCubicVolume = 1 ;

#if 0  
  // add test of dimensions 

   {
      G4cerr << "ERROR - G4TwistedGenTrap()::G4TwistedGenTrap(): " << GetName() << G4endl
             << "        Dimensions too small ! - "
             << pDx1 << ", " << pDx2 << ", " << pDy << ", " << pDz << G4endl 
             << " twistangle " << twistedangle/deg << " deg" << G4endl ;
      
      G4Exception("G4TwistedGenTrap::G4TwistedGenTrap()", "InvalidSetup",
                  FatalException, "Invalid dimensions. Too small, or twist angle too big.");
  }
#endif


}


//=====================================================================
//* destructor --------------------------------------------------------

G4TwistedGenTrap::~G4TwistedGenTrap()
{

  if (fLowerEndcap) delete fLowerEndcap ;
  if (fUpperEndcap) delete fUpperEndcap ;

  if (fSide0)      delete  fSide0 ;
  if (fSide90)     delete  fSide90 ;
  if (fSide180)    delete  fSide180 ;
  if (fSide270)    delete  fSide270 ;


}

//=====================================================================
//* ComputeDimensions -------------------------------------------------

void G4TwistedGenTrap::ComputeDimensions(G4VPVParameterisation* ,
                              const G4int ,
                              const G4VPhysicalVolume* )
{
  G4Exception("G4TwistedGenTrap::ComputeDimensions()",
              "NotSupported", FatalException,
              "G4TwistedGenTrap does not support Parameterisation.");

}

//=====================================================================
//* CalculateExtent ---------------------------------------------------

G4bool G4TwistedGenTrap::CalculateExtent( const EAxis        pAxis,
                                       const G4VoxelLimits     &pVoxelLimit,
                                       const G4AffineTransform &pTransform,
                                             G4double          &pMin,
                                             G4double          &pMax ) const
{

  G4double fDy = fDy1 ;  // temp !!!!!!!!!!!!

  G4double maxSide = ( fDx1 > fDx2 ? fDx1 : fDx2 ) ;
  G4double maxRad = std::sqrt( maxSide*maxSide + fDy*fDy);

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
               xMax < pVoxelLimit.GetMinXExtent()-kCarTolerance  ) return false ;
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
               yMax < pVoxelLimit.GetMinYExtent()-kCarTolerance   ) return false ;
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
               zMax < pVoxelLimit.GetMinZExtent()-kCarTolerance   ) return false ;
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
       ( pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5);
          
      if ( pMin != kInfinity || pMax != -kInfinity )
        {
          existsAfterClip = true ;
          

          // Check to see if endpoints are in the solid
          
          clipCentre(pAxis) = pVoxelLimit.GetMinExtent(pAxis);
          
          if (Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside)
            {
              pMin = pVoxelLimit.GetMinExtent(pAxis);
            }
          else
            {
              pMin -= kCarTolerance;
            }
          clipCentre(pAxis) = pVoxelLimit.GetMaxExtent(pAxis);
          
          if (Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside)
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

G4ThreeVectorList*
G4TwistedGenTrap::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{

  G4double fDy = fDy1 ;  // temp !!!!!!!!!!


  G4ThreeVectorList* vertices = new G4ThreeVectorList();
  vertices->reserve(8);

  if (vertices)
  {

    G4double maxSide = ( fDx1 > fDx2 ? fDx1 : fDx2 ) ;
    G4double maxRad = std::sqrt( maxSide*maxSide + fDy*fDy);

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
    G4Exception("G4TwistedGenTrap::CreateRotatedVertices()",
                "FatalError", FatalException,
                "Error in allocation of vertices. Out of memory !");
  }
  return vertices;
}

//=====================================================================
//* Inside ------------------------------------------------------------

EInside G4TwistedGenTrap::Inside(const G4ThreeVector& p) const
{

  G4double fDy = fDy1 ;  // temp !!!!!!!!!!!!!!!!!1

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

   G4double phi = - p.z()/(2*fDz) * fPhiTwist ;  // rotate the point to z=0
   G4double cphi = std::cos(phi) ;
   G4double sphi = std::sin(phi) ;
   G4double posx = p.x() * cphi - p.y() * sphi   ;
   G4double posy = p.x() * sphi + p.y() * cphi   ;
   G4double posz = p.z()  ;

   G4double wpos = fDx2 + ( fDx1-fDx2)/2. - posy * (fDx1-fDx2)/(2*fDy) ;


#ifdef G4SPECSDEBUG

   G4cout << "inside called: p = " << p << G4endl ; 
   G4cout << "Trapezoid is 1/2*(x,y,z) : " << fDx1 << ", " << fDx2 << ", " << 
   fDy << ", " << fDz << G4endl ;
   G4cout << "posx = " << posx << G4endl ;
   G4cout << "posy = " << posy << G4endl ;

  G4cout << "wpos = " << wpos << G4endl ;

#endif 

  if ( std::fabs(posx) <= wpos - kCarTolerance*0.5 )
  {
    if (std::fabs(posy) <= fDy - kCarTolerance*0.5 )
    {
      if      (std::fabs(posz) <= fDz - kCarTolerance*0.5 ) *tmpin = kInside ;
      else if (std::fabs(posz) <= fDz + kCarTolerance*0.5 ) *tmpin = kSurface ;
    }
    else if (std::fabs(posy) <= fDy + kCarTolerance*0.5 )
    {
      if (std::fabs(posz) <= fDz + kCarTolerance*0.5 ) *tmpin = kSurface ;
    }
  }
  else if (std::fabs(posx) <= wpos + kCarTolerance*0.5 )
  {
    if (std::fabs(posy) <= fDy + kCarTolerance*0.5 )
    {
      if (std::fabs(posz) <= fDz + kCarTolerance*0.5) *tmpin = kSurface ;
    }
  }
  
   return fLastInside.inside;

}

//=====================================================================
//* SurfaceNormal -----------------------------------------------------

G4ThreeVector G4TwistedGenTrap::SurfaceNormal(const G4ThreeVector& p) const
{
   //
   // return the normal unit vector to the Hyperbolical Surface at a point 
   // p on (or nearly on) the surface
   //
   // Which of the three or four surfaces are we closest to?
   //


   if (fLastNormal.p == p) {

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
   for (i=0; i< 6; i++) {
      G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
      if (tmpdistance < distance) {
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

G4double G4TwistedGenTrap::DistanceToIn (const G4ThreeVector& p,
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
   if (fLastDistanceToInWithV.p == p && fLastDistanceToInWithV.vec == v) {


     return fLastDistanceToIn.value;
   } else {
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

   if (currentside == kInside) {

   } else if (currentside == kSurface) {
      // particle is just on a boundary.
      // if the particle is entering to the volume, return 0.
      G4ThreeVector normal = SurfaceNormal(p);
      if (normal*v < 0) {

         *tmpdist = 0;
         return fLastDistanceToInWithV.value;
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
   
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   G4int i;
   G4int besti = -1;
   for (i=0; i < 6 ; i++) {

#ifdef G4SPECSDEBUG
      G4cout << G4endl << "surface " << i << ": " << G4endl << G4endl ;
#endif
      G4double tmpdistance = surfaces[i]->DistanceToIn(p, v, xx);
#ifdef G4SPECSDEBUG
      G4cout << "Solid DistanceToIn : distance = " << tmpdistance << G4endl ; 
      G4cout << "intersection point = " << xx << G4endl ;
#endif 
      if (tmpdistance < distance) {
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

#ifdef DISTANCETOIN
//=====================================================================
//* DistanceToIn (p) --------------------------------------------------

G4double G4TwistedGenTrap::DistanceToIn (const G4ThreeVector& p) const
{
   // DistanceToIn(p):
   // Calculate distance to surface of shape from `outside',
   // allowing for tolerance
   //

   //
   // checking last value
   //
   

  G4double fDy = fDy1 ;  // temp !!!!!!!!!!!!!!!!!1

   G4ThreeVector *tmpp;
   G4double      *tmpdist;
   if (fLastDistanceToIn.p == p) {
     return fLastDistanceToIn.value;
   } else {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToIn.p));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToIn.value));
      tmpp->set(p.x(), p.y(), p.z());
   }


   //
   // Calculate DistanceToIn(p) 
   //
   
   EInside currentside = Inside(p);

   switch (currentside) {

      case (kInside) : {
      }

      case (kSurface) : {
         *tmpdist = 0.;
         return fLastDistanceToIn.value;
      }

      case (kOutside) : {
         // Initialize

        G4double safex, safey, safez, safe = 0.0 ;

        G4double maxSide = ( fDx1 > fDx2 ? fDx1 : fDx2 ) ;
        G4double maxRad = std::sqrt( maxSide*maxSide + fDy*fDy);

        G4cout << "maxRad = " << maxRad << G4endl ;

        safex = std::fabs(p.x()) - maxRad ;
        safey = std::fabs(p.y()) - maxRad ;
        safez = std::fabs(p.z()) - fDz ;
        
        if (safex > safe) safe = safex ;
        if (safey > safe) safe = safey ;
        if (safez > safe) safe = safez ;
        
        *tmpdist = safe;        
        return fLastDistanceToIn.value;
        
      }


      default : {
         G4Exception("G4TwistedGenTrap::DistanceToIn(p)", "InvalidCondition",
                     FatalException, "Unknown point location!");
      }
   } // switch end
   return kInfinity;
}
#else

G4double G4TwistedGenTrap::DistanceToIn (const G4ThreeVector& p) const
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
   if (fLastDistanceToIn.p == p) {

     return fLastDistanceToIn.value;
   } else {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToIn.p));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToIn.value));
      tmpp->set(p.x(), p.y(), p.z());
   }

   //
   // Calculate DistanceToIn(p) 
   //
   
   EInside currentside = Inside(p);

   switch (currentside) {

      case (kInside) : {

      }

      case (kSurface) : {
         *tmpdist = 0.;
         return fLastDistanceToIn.value;
      }

      case (kOutside) : {
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
         for (i=0; i< 6; i++) {
            G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
            if (tmpdistance < distance)  {
               distance = tmpdistance;
               bestxx = xx;
               besti = i;
            }
         }
      
         *tmpdist = distance;
         return fLastDistanceToIn.value;
      }

      default : {
         G4Exception("G4TwistedGenTrap::DistanceToIn(p)", "InvalidCondition",
                     FatalException, "Unknown point location!");
      }
   } // switch end
   return kInfinity;
}

#endif

//=====================================================================
//* DistanceToOut (p, v) ----------------------------------------------

G4double G4TwistedGenTrap::DistanceToOut( const G4ThreeVector& p,
                                       const G4ThreeVector& v,
                                       const G4bool calcNorm,
                                       G4bool *validNorm, G4ThreeVector *norm ) const
{
   // DistanceToOut (p, v):
   // Calculate distance to surface of shape from `inside'
   // along with the v, allowing for tolerance.
   // The function returns kInfinity if no intersection or
   // just grazing within tolerance.
   //
   

   //
   // checking last value
   //
   
   G4ThreeVector *tmpp;
   G4ThreeVector *tmpv;
   G4double      *tmpdist;
   if (fLastDistanceToOutWithV.p == p && fLastDistanceToOutWithV.vec == v  ) {
      return fLastDistanceToOutWithV.value;
   } else {
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

   if (currentside == kOutside) {

   } else if (currentside == kSurface) {
      // particle is just on a boundary.
      // if the particle is exiting from the volume, return 0.
      G4ThreeVector normal = SurfaceNormal(p);
      G4VSurface *blockedsurface = fLastNormal.surface[0];
      if (normal*v > 0) {
            if (calcNorm) {
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
      if (tmpdistance < distance) {
         distance = tmpdistance;
         bestxx = xx; 
         besti = i;
      }
   }

   if (calcNorm) {
      if (besti != -1) {
         *norm = (surfaces[besti]->GetNormal(p, true));
         *validNorm = surfaces[besti]->IsValidNorm();
      }
   }

   *tmpdist = distance;
   // timer.Stop();
   return fLastDistanceToOutWithV.value;
}

#ifdef DISTANCETOIN
//=====================================================================
//* DistanceToOut (p) ----------------------------------------------

G4double G4TwistedGenTrap::DistanceToOut( const G4ThreeVector& p ) const
{
   // DistanceToOut(p):
   // Calculate distance to surface of shape from `inside', 
   // allowing for tolerance
   //
   
   //
   // checking last value
   //

  G4double fDy = fDy1 ;

   G4ThreeVector *tmpp;
   G4double      *tmpdist;
   if (fLastDistanceToOut.p == p) {
      return fLastDistanceToOut.value;
   } else {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToOut.p));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToOut.value));
      tmpp->set(p.x(), p.y(), p.z());
   }
   
   //
   // Calculate DistanceToOut(p)
   //
   
   EInside currentside = Inside(p);

   switch (currentside) {
      case (kOutside) : {

      }
      case (kSurface) : {
        *tmpdist = 0.;
         return fLastDistanceToOut.value;
      }
      
      case (kInside) : {
         // Initialize

        G4double      rminX = ( fDx1 < fDx2  ? fDx1 : fDx2 )  ;
        G4double      rmin = (  rminX < fDy  ? rminX : fDy )  ;

        G4double safx1,safx2,safy1,safy2,safz1,safz2,safe=0.0;

        safx1 = rmin - p.x() ;
        safx2 = rmin + p.x() ;
        safy1 = rmin - p.y() ;
        safy2 = rmin + p.y() ;
        safz1 = fDz - p.z() ;
        safz2 = fDz + p.z() ;  
        
  // shortest Dist to any boundary now MIN(safx1,safx2,safy1..)

        if (safx2 < safx1) safe = safx2 ;
        else               safe = safx1 ;
        if (safy1 < safe)  safe = safy1 ;
        if (safy2 < safe)  safe = safy2 ;
        if (safz1 < safe)  safe = safz1 ;
        if (safz2 < safe)  safe = safz2 ;
        
        if (safe < 0) safe = 0 ;

        *tmpdist = safe;

        return fLastDistanceToOut.value;
      }
      
      default : {
        G4Exception("G4TwistedGenTrap::DistanceToOut(p)", "InvalidCondition",
                 FatalException, "Unknown point location!");
      }

   } // switch end

   return 0;

}

#else

G4double G4TwistedGenTrap::DistanceToOut( const G4ThreeVector& p ) const
{
   // DistanceToOut(p):
   // Calculate distance to surface of shape from `inside', 
   // allowing for tolerance
   //


   //
   // checking last value
   //

   
   G4ThreeVector *tmpp;
   G4double      *tmpdist;
   if (fLastDistanceToOut.p == p) {
      return fLastDistanceToOut.value;
   } else {
      tmpp    = const_cast<G4ThreeVector*>(&(fLastDistanceToOut.p));
      tmpdist = const_cast<G4double*>(&(fLastDistanceToOut.value));
      tmpp->set(p.x(), p.y(), p.z());
   }
   
   //
   // Calculate DistanceToOut(p)
   //
   
   EInside currentside = Inside(p);

   switch (currentside) {
      case (kOutside) : {

      }
      case (kSurface) : {
        *tmpdist = 0.;
         return fLastDistanceToOut.value;
      }
      
      case (kInside) : {
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
         for (i=0; i< 6; i++) {
            G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
            if (tmpdistance < distance) {
               distance = tmpdistance;
               bestxx = xx;
               besti = i;
            }
         }


         *tmpdist = distance;
   
         return fLastDistanceToOut.value;
      }
      
      default : {
         G4Exception("G4TwistedGenTrap::DistanceToOut(p)", "InvalidCondition",
                     FatalException, "Unknown point location!");
      }
   } // switch end
   return 0;
}


#endif

//=====================================================================
//* StreamInfo --------------------------------------------------------

std::ostream& G4TwistedGenTrap::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4TwistedGenTrap\n"
     << " Parameters: \n"
     << "    Twist angle   : " << fPhiTwist/degree << " degrees \n"
     << "-----------------------------------------------------------\n";

  return os;
}


//=====================================================================
//* DiscribeYourselfTo ------------------------------------------------

void G4TwistedGenTrap::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddThis (*this);
}

//=====================================================================
//* GetExtent ---------------------------------------------------------

G4VisExtent G4TwistedGenTrap::GetExtent() const 
{
  G4double fDy = fDy1 ;  // temp !!!!!!!!!!!!!!!!!!!!!1

  G4double maxSide = ( fDx1 > fDx2 ? fDx1 : fDx2 ) ;
  G4double maxRad = std::sqrt( maxSide*maxSide + fDy*fDy);

  return G4VisExtent(-maxRad, maxRad ,
                     -maxRad, maxRad ,
                     -fDz, fDz );
}

//=====================================================================
//* CreatePolyhedron --------------------------------------------------

G4Polyhedron* G4TwistedGenTrap::CreatePolyhedron () const 
{
   //  return new G4PolyhedronHype (fRMin, fRMax, fDz, fSPhi, fDPhi);
   return 0;
}

//=====================================================================
//* CreateNUBS --------------------------------------------------------

G4NURBS* G4TwistedGenTrap::CreateNURBS () const 
{

  G4double fDy = fDy1 ;  // temp !!!!!!!!!!!!!!!!!!!!!1

  G4double maxSide = ( fDx1 > fDx2 ? fDx1 : fDx2 ) ;
  G4double maxRad = std::sqrt( maxSide*maxSide + fDy*fDy);

  return new G4NURBStube(maxRad, maxRad, fDz); 
   // Tube for now!!!
}

//=====================================================================
//* CreateSurfaces ----------------------------------------------------

void G4TwistedGenTrap::CreateSurfaces() 
{
   
   // create 6 surfaces of TwistedTub.

  G4double fDy = fDy1 ;  // temp !!!!!!!!!!!!!!!!!!!!!!1

   fSide0   = new G4TwistedTrapAlphaSide("0deg"   ,fPhiTwist, fDz, fTheta,  fPhi    , fDy1, fDx1, fDx2, fDy2, fDx3, fDx4, fAlph, 0.*deg ) ;
   fSide180 = new G4TwistedTrapAlphaSide("180deg", fPhiTwist, fDz, fTheta, fPhi + pi, fDy1, fDx2, fDx1, fDy2, fDx4, fDx3, fAlph, 180.*deg) ;

   fSide180 = new G4TwistedTrapParallelSide("90deg",  fPhiTwist, fDz, fTheta, fPhi,      fDy1, fDx1, fDx2, fDy2, fDx3, fDx4, fAlph, 0.*deg) ;
   fSide180 = new G4TwistedTrapParallelSide("270deg", fPhiTwist, fDz, fTheta, fPhi + pi, fDy1, fDx2, fDx1, fDy2, fDx4, fDx3, fAlph, 180.*deg) ;


   fUpperEndcap = new G4FlatTrapSide("UpperCap",fPhiTwist, fDx1, fDx2, fDy, fDz,  1 ) ;
   fLowerEndcap = new G4FlatTrapSide("LowerCap",fPhiTwist, fDx1, fDx2, fDy, fDz, -1 ) ;
 
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

G4GeometryType G4TwistedGenTrap::GetEntityType() const
{
  return G4String("G4TwistedGenTrap");
}
