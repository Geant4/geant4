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
//
// $Id: G4TwistedTrapezoid.cc,v 1.2 2004-08-27 13:36:51 link Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistedTrapezoid.cc
//
// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------

#include "G4TwistedTrapezoid.hh"

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

G4TwistedTrapezoid::G4TwistedTrapezoid(const G4String &pname,
				       G4double  twistedangle,
				       G4double  halfzlen,
				       G4double  halfSideX,
				       G4double  halfSideY,
				       G4double  dphi)
  : G4VSolid(pname), 
     fLowerEndcap(0), fUpperEndcap(0), fSide0(0),
    fSide90(0), fSide180(0), fSide270(0)
{

   SetFields(twistedangle, halfzlen, halfSideX, halfSideY);
   CreateSurfaces();

}



//=====================================================================
//* destructor --------------------------------------------------------

G4TwistedTrapezoid::~G4TwistedTrapezoid()
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

void G4TwistedTrapezoid::ComputeDimensions(G4VPVParameterisation* /* p */ ,
				      const G4int            /* n  */ ,
				      const G4VPhysicalVolume* /* pRep */ )
{
  G4Exception("G4TwistedTrapezoid::ComputeDimensions()",
              "NotSupported", FatalException,
	      "G4TwistedTrapezoid does not support Parameterisation.");
}


//=====================================================================
//* CalculateExtent ---------------------------------------------------

G4bool G4TwistedTrapezoid::CalculateExtent( const EAxis              axis,
                                       const G4VoxelLimits     &voxelLimit,
                                       const G4AffineTransform &transform,
                                             G4double          &min,
                                             G4double          &max ) const
{

  // missing 
  G4cout << "Calculate Extent missing\n" << G4endl ;
  return true ;

}


//=====================================================================
//* AddPolyToExtent ---------------------------------------------------

void G4TwistedTrapezoid::AddPolyToExtent( const G4ThreeVector &v0,
                                     const G4ThreeVector &v1,
                                     const G4ThreeVector &w1,
                                     const G4ThreeVector &w0,
                                     const G4VoxelLimits &voxelLimit,
                                     const EAxis          axis,
                                     G4SolidExtentList   &extentList ) 
{
    // Utility function for CalculateExtent
    //
    G4ClippablePolygon phiPoly;

    phiPoly.AddVertexInOrder( v0 );
    phiPoly.AddVertexInOrder( v1 );
    phiPoly.AddVertexInOrder( w1 );
    phiPoly.AddVertexInOrder( w0 );

    if (phiPoly.PartialClip( voxelLimit, axis )) {
        phiPoly.SetNormal( (v1-v0).cross(w0-v0).unit() );
        extentList.AddSurface( phiPoly );
    }
}


//=====================================================================
//* Inside ------------------------------------------------------------

EInside G4TwistedTrapezoid::Inside(const G4ThreeVector& p) const
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

   G4double phi = - p.z()/(2*fZHalfLength) * fPhiTwist ;
   G4double cphi = cos(phi) ;
   G4double sphi = sin(phi) ;
   G4double posx = p.x() * cphi - p.y() * sphi   ;
   G4double posy = p.x() * sphi + p.y() * cphi   ;
   G4double posz = p.z()  ;

   G4cout << "inside called: p = " << p << G4endl ; 
   G4cout << "Trapezoid is : " << fHalfSides[0] << ", " << 
   fHalfSides[1] << ", " << fZHalfLength << G4endl ;

   G4cout << "distanceToOut = ( " << 
    posx << ", " <<
    posy << ", " <<
    posz <<") " <<
    G4endl ;

  if ( fabs(posx) <= fHalfSides[0] - kCarTolerance*0.5 )
  {
    if (fabs(posy) <= fHalfSides[1] - kCarTolerance*0.5 )
    {
      if      (fabs(posz) <= fZHalfLength - kCarTolerance*0.5 ) *tmpin = kInside ;
      else if (fabs(posz) <= fZHalfLength + kCarTolerance*0.5 ) *tmpin = kSurface ;
    }
    else if (fabs(posy) <= fHalfSides[1] + kCarTolerance*0.5 )
    {
      if (fabs(posz) <= fZHalfLength + kCarTolerance*0.5 ) *tmpin = kSurface ;
    }
  }
  else if (fabs(posx) <= fHalfSides[0] + kCarTolerance*0.5 )
  {
    if (fabs(posy) <= fHalfSides[1] + kCarTolerance*0.5 )
    {
      if (fabs(posz) <= fZHalfLength + kCarTolerance*0.5) *tmpin = kSurface ;
    }
  }
  
   return fLastInside.inside;

}

//=====================================================================
//* SurfaceNormal -----------------------------------------------------

G4ThreeVector G4TwistedTrapezoid::SurfaceNormal(const G4ThreeVector& p) const
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

   /*
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;
   */

   G4ThreeVector xx;
   G4ThreeVector bestxx;
   G4int i;
   G4int besti = -1;
   for (i=0; i< 4; i++) {
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

G4double G4TwistedTrapezoid::DistanceToIn (const G4ThreeVector& p,
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

   /*
   surfaces[0] = fLowerEndcap;
   surfaces[1] = fUpperEndcap;
   surfaces[2] = fLatterTwisted;
   surfaces[3] = fFormerTwisted;
   */
   //   surfaces[4] = fInnerHype;
   //surfaces[5] = fOuterHype;
   
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   G4int i;
   G4int besti = -1;
   for (i=0; i< 4; i++) {
      G4double tmpdistance = surfaces[i]->DistanceToIn(p, v, xx);
      if (tmpdistance < distance) {
         distance = tmpdistance;
         bestxx = xx;
         besti = i;
      }
   }

   *tmpdist = distance;
   // timer.Stop();
   return fLastDistanceToInWithV.value;
}
 
//=====================================================================
//* DistanceToIn (p) --------------------------------------------------

G4double G4TwistedTrapezoid::DistanceToIn (const G4ThreeVector& p) const
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

	 /*
         surfaces[4] = fLatterTwisted;
         surfaces[5] = fFormerTwisted;
	 */

         G4int i;
         G4int besti = -1;
         G4ThreeVector xx;
         G4ThreeVector bestxx;
         for (i=0; i< 4; i++) {
            G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
            if (tmpdistance < distance)  {
               distance = tmpdistance;
               bestxx = xx;
               besti = i;
            }
         }
      
         *tmpdist = distance;
         // timer.Stop();
         return fLastDistanceToIn.value;
      }

      default : {
         G4Exception("G4TwistedTrapezoid::DistanceToIn(p)", "InvalidCondition",
                     FatalException, "Unknown point location!");
      }
   } // switch end
   return kInfinity;
}

//=====================================================================
//* DistanceToOut (p, v) ----------------------------------------------

G4double G4TwistedTrapezoid::DistanceToOut( const G4ThreeVector& p,
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

   /*
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;
   */
   G4int i;
   G4int besti = -1;
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   for (i=0; i< 4; i++) {
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


//=====================================================================
//* DistanceToOut (p) ----------------------------------------------

G4double G4TwistedTrapezoid::DistanceToOut( const G4ThreeVector& p ) const
{
   // DistanceToOut(p):
   // Calculate distance to surface of shape from `inside', 
   // allowing for tolerance
   //
   
   //
   // opening message
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

	 /*
         surfaces[4] = fLowerEndcap;
         surfaces[5] = fUpperEndcap;
	 */
         G4int i;
         G4int besti = -1;
         G4ThreeVector xx;
         G4ThreeVector bestxx;
         for (i=0; i< 4; i++) {
            G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
            if (tmpdistance < distance) {
               distance = tmpdistance;
               bestxx = xx;
               besti = i;
            }
         }

         *tmpdist = distance;
         // timer.Stop();
         return fLastDistanceToOut.value;
      }
      
      default : {
         G4Exception("G4TwistedTrapezoid::DistanceToOut(p)", "InvalidCondition",
                     FatalException, "Unknown point location!");
      }
   } // switch end
   return 0;
}

//=====================================================================
//* StreamInfo --------------------------------------------------------

std::ostream& G4TwistedTrapezoid::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4TwistedTrapezoid\n"
     << " Parameters: \n"
     << "    Twist angle   : " << fPhiTwist/degree << " degrees \n"
     << "-----------------------------------------------------------\n";

  return os;
}


//=====================================================================
//* DiscribeYourselfTo ------------------------------------------------

void G4TwistedTrapezoid::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddThis (*this);
}

//=====================================================================
//* GetExtent ---------------------------------------------------------

G4VisExtent G4TwistedTrapezoid::GetExtent() const 
{

  return G4VisExtent(-fHalfSides[0], fHalfSides[0],
		     -fHalfSides[1], fHalfSides[1],
		     -fZHalfLength, fZHalfLength );
}

//=====================================================================
//* CreatePolyhedron --------------------------------------------------

G4Polyhedron* G4TwistedTrapezoid::CreatePolyhedron () const 
{
   //  return new G4PolyhedronHype (fRMin, fRMax, fDz, fSPhi, fDPhi);
   return 0;
}

//=====================================================================
//* CreateNUBS --------------------------------------------------------

G4NURBS* G4TwistedTrapezoid::CreateNURBS () const 
{
  G4double maxRad = sqrt( fHalfSides[0]*fHalfSides[0] + fHalfSides[1]*fHalfSides[1]);
   return new G4NURBStube(maxRad, maxRad, fZHalfLength); 
   // Tube for now!!!
}

//=====================================================================
//* CreateSurfaces ----------------------------------------------------

void G4TwistedTrapezoid::CreateSurfaces() 
{
   
   // create 6 surfaces of TwistedTub.

   fSide0 = new G4TwistedTrapSide("0deg", fPhiTwist, fZHalfLength, fHalfSides, 0) ;
   fSide90 = new G4TwistedTrapSide("90deg", fPhiTwist, fZHalfLength, fHalfSides, 90) ;
   fSide180 = new G4TwistedTrapSide("180deg", fPhiTwist, fZHalfLength, fHalfSides, 180) ;
   fSide270 = new G4TwistedTrapSide("270deg", fPhiTwist, fZHalfLength, fHalfSides, 270) ;

}


//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4TwistedTrapezoid::GetEntityType() const
{
  return G4String("G4TwistedTrapezoid");
}
