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
// G4TwistedTubs implementation
//
// 01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp), created.
// 13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//               from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------

#include "G4TwistedTubs.hh"

#include "G4GeomTools.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"
#include "G4ClippablePolygon.hh"
#include "G4VPVParameterisation.hh"
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

G4TwistedTubs::G4TwistedTubs(const G4String& pname,
                                   G4double  twistedangle,
                                   G4double  endinnerrad,
                                   G4double  endouterrad,
                                   G4double  halfzlen,
                                   G4double  dphi)
   : G4VSolid(pname), fDPhi(dphi), 
     fLowerEndcap(0), fUpperEndcap(0), fLatterTwisted(0),
     fFormerTwisted(0), fInnerHype(0), fOuterHype(0)
{
   if (endinnerrad < DBL_MIN)
   {
      G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002",
                  FatalErrorInArgument, "Invalid end-inner-radius!");
   }
            
   G4double sinhalftwist = std::sin(0.5 * twistedangle);

   G4double endinnerradX = endinnerrad * sinhalftwist;
   G4double innerrad     = std::sqrt( endinnerrad * endinnerrad
                                 - endinnerradX * endinnerradX );

   G4double endouterradX = endouterrad * sinhalftwist;
   G4double outerrad     = std::sqrt( endouterrad * endouterrad
                                 - endouterradX * endouterradX );
   
   // temporary treatment!!
   SetFields(twistedangle, innerrad, outerrad, -halfzlen, halfzlen);
   CreateSurfaces();
}

G4TwistedTubs::G4TwistedTubs(const G4String& pname,     
                                   G4double  twistedangle,    
                                   G4double  endinnerrad,  
                                   G4double  endouterrad,  
                                   G4double  halfzlen,
                                   G4int     nseg,
                                   G4double  totphi)
   : G4VSolid(pname),
     fLowerEndcap(0), fUpperEndcap(0), fLatterTwisted(0),
     fFormerTwisted(0), fInnerHype(0), fOuterHype(0)
{

   if (!nseg)
   {
      std::ostringstream message;
      message << "Invalid number of segments." << G4endl
              << "        nseg = " << nseg;
      G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002",
                  FatalErrorInArgument, message);
   }
   if (totphi == DBL_MIN || endinnerrad < DBL_MIN)
   {
      G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002",
                FatalErrorInArgument, "Invalid total-phi or end-inner-radius!");
   }
         
   G4double sinhalftwist = std::sin(0.5 * twistedangle);

   G4double endinnerradX = endinnerrad * sinhalftwist;
   G4double innerrad     = std::sqrt( endinnerrad * endinnerrad
                                 - endinnerradX * endinnerradX );

   G4double endouterradX = endouterrad * sinhalftwist;
   G4double outerrad     = std::sqrt( endouterrad * endouterrad
                                 - endouterradX * endouterradX );
   
   // temporary treatment!!
   fDPhi = totphi / nseg;
   SetFields(twistedangle, innerrad, outerrad, -halfzlen, halfzlen);
   CreateSurfaces();
}

G4TwistedTubs::G4TwistedTubs(const G4String& pname,
                                   G4double  twistedangle,
                                   G4double  innerrad,
                                   G4double  outerrad,
                                   G4double  negativeEndz,
                                   G4double  positiveEndz,
                                   G4double  dphi)
   : G4VSolid(pname), fDPhi(dphi),
     fLowerEndcap(0), fUpperEndcap(0), fLatterTwisted(0),
     fFormerTwisted(0), fInnerHype(0), fOuterHype(0)
{
   if (innerrad < DBL_MIN)
   {
      G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002",
                  FatalErrorInArgument, "Invalid end-inner-radius!");
   }
                 
   SetFields(twistedangle, innerrad, outerrad, negativeEndz, positiveEndz);
   CreateSurfaces();
}

G4TwistedTubs::G4TwistedTubs(const G4String& pname,     
                                   G4double  twistedangle,    
                                   G4double  innerrad,  
                                   G4double  outerrad,  
                                   G4double  negativeEndz,
                                   G4double  positiveEndz,
                                   G4int     nseg,
                                   G4double  totphi)
   : G4VSolid(pname),
     fLowerEndcap(0), fUpperEndcap(0), fLatterTwisted(0),
     fFormerTwisted(0), fInnerHype(0), fOuterHype(0)
{
   if (!nseg)
   {
      std::ostringstream message;
      message << "Invalid number of segments." << G4endl
              << "        nseg = " << nseg;
      G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002",
                  FatalErrorInArgument, message);
   }
   if (totphi == DBL_MIN || innerrad < DBL_MIN)
   {
      G4Exception("G4TwistedTubs::G4TwistedTubs()", "GeomSolids0002",
                FatalErrorInArgument, "Invalid total-phi or end-inner-radius!");
   }
            
   fDPhi = totphi / nseg;
   SetFields(twistedangle, innerrad, outerrad, negativeEndz, positiveEndz);
   CreateSurfaces();
}

//=====================================================================
//* Fake default constructor ------------------------------------------

G4TwistedTubs::G4TwistedTubs( __void__& a )
  : G4VSolid(a), fPhiTwist(0.), fInnerRadius(0.), fOuterRadius(0.), fDPhi(0.),
    fZHalfLength(0.), fInnerStereo(0.), fOuterStereo(0.), fTanInnerStereo(0.),
    fTanOuterStereo(0.), fKappa(0.), fInnerRadius2(0.), fOuterRadius2(0.),
    fTanInnerStereo2(0.), fTanOuterStereo2(0.), fLowerEndcap(0), fUpperEndcap(0),
    fLatterTwisted(0), fFormerTwisted(0), fInnerHype(0), fOuterHype(0)
{
}

//=====================================================================
//* destructor --------------------------------------------------------

G4TwistedTubs::~G4TwistedTubs()
{
   if (fLowerEndcap)   { delete fLowerEndcap;   }
   if (fUpperEndcap)   { delete fUpperEndcap;   }
   if (fLatterTwisted) { delete fLatterTwisted; }
   if (fFormerTwisted) { delete fFormerTwisted; }
   if (fInnerHype)     { delete fInnerHype;     }
   if (fOuterHype)     { delete fOuterHype;     }
   if (fpPolyhedron)   { delete fpPolyhedron; fpPolyhedron = nullptr; }
}

//=====================================================================
//* Copy constructor --------------------------------------------------

G4TwistedTubs::G4TwistedTubs(const G4TwistedTubs& rhs)
  : G4VSolid(rhs), fPhiTwist(rhs.fPhiTwist),
    fInnerRadius(rhs.fInnerRadius), fOuterRadius(rhs.fOuterRadius),
    fDPhi(rhs.fDPhi), fZHalfLength(rhs.fZHalfLength),
    fInnerStereo(rhs.fInnerStereo), fOuterStereo(rhs.fOuterStereo),
    fTanInnerStereo(rhs.fTanInnerStereo), fTanOuterStereo(rhs.fTanOuterStereo),
    fKappa(rhs.fKappa), fInnerRadius2(rhs.fInnerRadius2), 
    fOuterRadius2(rhs.fOuterRadius2), fTanInnerStereo2(rhs.fTanInnerStereo2),
    fTanOuterStereo2(rhs.fTanOuterStereo2),
    fLowerEndcap(0), fUpperEndcap(0), fLatterTwisted(0), fFormerTwisted(0),
    fInnerHype(0), fOuterHype(0),
    fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
    fLastInside(rhs.fLastInside), fLastNormal(rhs.fLastNormal),
    fLastDistanceToIn(rhs.fLastDistanceToIn),
    fLastDistanceToOut(rhs.fLastDistanceToOut),
    fLastDistanceToInWithV(rhs.fLastDistanceToInWithV),
    fLastDistanceToOutWithV(rhs.fLastDistanceToOutWithV)
{
  for (auto i=0; i<2; ++i)
  {
    fEndZ[i] = rhs.fEndZ[i];
    fEndInnerRadius[i] = rhs.fEndInnerRadius[i];
    fEndOuterRadius[i] = rhs.fEndOuterRadius[i];
    fEndPhi[i] = rhs.fEndPhi[i];
    fEndZ2[i] = rhs.fEndZ2[i];
  }
  CreateSurfaces();
}


//=====================================================================
//* Assignment operator -----------------------------------------------

G4TwistedTubs& G4TwistedTubs::operator = (const G4TwistedTubs& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   fPhiTwist= rhs.fPhiTwist;
   fInnerRadius= rhs.fInnerRadius; fOuterRadius= rhs.fOuterRadius;
   fDPhi= rhs.fDPhi; fZHalfLength= rhs.fZHalfLength;
   fInnerStereo= rhs.fInnerStereo; fOuterStereo= rhs.fOuterStereo;
   fTanInnerStereo= rhs.fTanInnerStereo; fTanOuterStereo= rhs.fTanOuterStereo;
   fKappa= rhs.fKappa; fInnerRadius2= rhs.fInnerRadius2; 
   fOuterRadius2= rhs.fOuterRadius2; fTanInnerStereo2= rhs.fTanInnerStereo2;
   fTanOuterStereo2= rhs.fTanOuterStereo2;
   fLowerEndcap= fUpperEndcap= fLatterTwisted= fFormerTwisted= 0;
   fInnerHype= fOuterHype= 0;
   fCubicVolume= rhs.fCubicVolume; fSurfaceArea= rhs.fSurfaceArea;
   fLastInside= rhs.fLastInside; fLastNormal= rhs.fLastNormal;
   fLastDistanceToIn= rhs.fLastDistanceToIn;
   fLastDistanceToOut= rhs.fLastDistanceToOut;
   fLastDistanceToInWithV= rhs.fLastDistanceToInWithV;
   fLastDistanceToOutWithV= rhs.fLastDistanceToOutWithV;
 
   for (auto i=0; i<2; ++i)
   {
     fEndZ[i] = rhs.fEndZ[i];
     fEndInnerRadius[i] = rhs.fEndInnerRadius[i];
     fEndOuterRadius[i] = rhs.fEndOuterRadius[i];
     fEndPhi[i] = rhs.fEndPhi[i];
     fEndZ2[i] = rhs.fEndZ2[i];
   }
 
   CreateSurfaces();
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = nullptr;

   return *this;
}

//=====================================================================
//* ComputeDimensions -------------------------------------------------

void G4TwistedTubs::ComputeDimensions(G4VPVParameterisation* /* p */ ,
                                      const G4int            /* n  */ ,
                                      const G4VPhysicalVolume* /* pRep */ )
{
  G4Exception("G4TwistedTubs::ComputeDimensions()",
              "GeomSolids0001", FatalException,
              "G4TwistedTubs does not support Parameterisation.");
}

//=====================================================================
//* BoundingLimits ----------------------------------------------------

void G4TwistedTubs::BoundingLimits(G4ThreeVector& pMin,
                                   G4ThreeVector& pMax) const
{
  // Find bounding tube
  G4double rmin = GetInnerRadius();
  G4double rmax = GetEndOuterRadius();

  G4double zmin = std::min(GetEndZ(0), GetEndZ(1));
  G4double zmax = std::max(GetEndZ(0), GetEndZ(1));

  G4double dphi = 0.5*GetDPhi();
  G4double sphi = std::min(GetEndPhi(0), GetEndPhi(1)) - dphi;
  G4double ephi = std::max(GetEndPhi(0), GetEndPhi(1)) + dphi;
  G4double totalphi = ephi - sphi;

  // Find bounding box
  if (dphi <= 0 || totalphi >= CLHEP::twopi)
  {
    pMin.set(-rmax,-rmax, zmin);
    pMax.set( rmax, rmax, zmax);
  }
  else
  {
    G4TwoVector vmin,vmax;
    G4GeomTools::DiskExtent(rmin, rmax, sphi, totalphi, vmin, vmax);
    pMin.set(vmin.x(), vmin.y(), zmin);
    pMax.set(vmax.x(), vmax.y(), zmax);
  }

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4TwistedTubs::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//=====================================================================
//* CalculateExtent ---------------------------------------------------

G4bool
G4TwistedTubs::CalculateExtent(const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                                     G4double& pMin, G4double& pMax) const
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

EInside G4TwistedTubs::Inside(const G4ThreeVector& p) const
{

   const G4double halftol
     = 0.5 * G4GeometryTolerance::GetInstance()->GetRadialTolerance();
   // static G4int timerid = -1;
   // G4Timer timer(timerid, "G4TwistedTubs", "Inside");
   // timer.Start();

   G4ThreeVector *tmpp;
   EInside       *tmpinside;
   if (fLastInside.p == p)
   {
     return fLastInside.inside;
   }
   else
   {
      tmpp      = const_cast<G4ThreeVector*>(&(fLastInside.p));
      tmpinside = const_cast<EInside*>(&(fLastInside.inside));
      tmpp->set(p.x(), p.y(), p.z());
   }
   
   EInside  outerhypearea = ((G4TwistTubsHypeSide *)fOuterHype)->Inside(p);
   G4double innerhyperho  = ((G4TwistTubsHypeSide *)fInnerHype)->GetRhoAtPZ(p);
   G4double distanceToOut = p.getRho() - innerhyperho; // +ve: inside

   if ((outerhypearea == kOutside) || (distanceToOut < -halftol))
   {
      *tmpinside = kOutside;
   }
   else if (outerhypearea == kSurface)
   {
      *tmpinside = kSurface;
   }
   else
   {
      if (distanceToOut <= halftol)
      {
         *tmpinside = kSurface;
      }
      else
      {
         *tmpinside = kInside;
      }
   }

   return fLastInside.inside;
}

//=====================================================================
//* SurfaceNormal -----------------------------------------------------

G4ThreeVector G4TwistedTubs::SurfaceNormal(const G4ThreeVector& p) const
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
   G4ThreeVector *tmpp          =
     const_cast<G4ThreeVector*>(&(fLastNormal.p));
   G4ThreeVector *tmpnormal     =
     const_cast<G4ThreeVector*>(&(fLastNormal.vec));
   G4VTwistSurface **tmpsurface =
     const_cast<G4VTwistSurface**>(fLastNormal.surface);
   tmpp->set(p.x(), p.y(), p.z());

   G4double      distance = kInfinity;

   G4VTwistSurface *surfaces[6];
   surfaces[0] = fLatterTwisted;
   surfaces[1] = fFormerTwisted;
   surfaces[2] = fInnerHype;
   surfaces[3] = fOuterHype;
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;

   G4ThreeVector xx;
   G4ThreeVector bestxx;
   G4int besti = -1;
   for (auto i=0; i<6; ++i)
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

G4double G4TwistedTubs::DistanceToIn (const G4ThreeVector& p,
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
   
   G4ThreeVector* tmpp;
   G4ThreeVector* tmpv;
   G4double* tmpdist;
   if ((fLastDistanceToInWithV.p == p) && (fLastDistanceToInWithV.vec == v))
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
   else
   {
     if (currentside == kSurface)
     {
       // particle is just on a boundary.
       // If the particle is entering to the volume, return 0.
       //
       G4ThreeVector normal = SurfaceNormal(p);
       if (normal*v < 0)
       {
         *tmpdist = 0.;
         return fLastDistanceToInWithV.value;
       } 
     }
   }
      
   // now, we can take smallest positive distance.
   
   // Initialize
   //
   G4double      distance = kInfinity;   

   // find intersections and choose nearest one.
   //
   G4VTwistSurface* surfaces[6];
   surfaces[0] = fLowerEndcap;
   surfaces[1] = fUpperEndcap;
   surfaces[2] = fLatterTwisted;
   surfaces[3] = fFormerTwisted;
   surfaces[4] = fInnerHype;
   surfaces[5] = fOuterHype;
   
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   for (auto i=0; i<6; ++i)
   {
      G4double tmpdistance = surfaces[i]->DistanceToIn(p, v, xx);
      if (tmpdistance < distance)
      {
         distance = tmpdistance;
         bestxx = xx;
      }
   }
   *tmpdist = distance;

   return fLastDistanceToInWithV.value;
}
 
//=====================================================================
//* DistanceToIn (p) --------------------------------------------------

G4double G4TwistedTubs::DistanceToIn (const G4ThreeVector& p) const
{
   // DistanceToIn(p):
   // Calculate distance to surface of shape from `outside',
   // allowing for tolerance
   
   //
   // checking last value
   //
   
   G4ThreeVector* tmpp;
   G4double* tmpdist;
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
      {}
      case (kSurface) :
      {
         *tmpdist = 0.;
         return fLastDistanceToIn.value;
      }
      case (kOutside) :
      {
         // Initialize
         G4double distance = kInfinity;   

         // find intersections and choose nearest one.
         G4VTwistSurface *surfaces[6];
         surfaces[0] = fLowerEndcap;
         surfaces[1] = fUpperEndcap;
         surfaces[2] = fLatterTwisted;
         surfaces[3] = fFormerTwisted;
         surfaces[4] = fInnerHype;
         surfaces[5] = fOuterHype;

         G4ThreeVector xx;
         G4ThreeVector bestxx;
         for (auto i=0; i<6; ++i)
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
         G4Exception("G4TwistedTubs::DistanceToIn(p)", "GeomSolids0003",
                     FatalException, "Unknown point location!");
      }
   } // switch end

   return kInfinity;
}

//=====================================================================
//* DistanceToOut (p, v) ----------------------------------------------

G4double G4TwistedTubs::DistanceToOut( const G4ThreeVector& p,
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
   
   G4ThreeVector* tmpp;
   G4ThreeVector* tmpv;
   G4double* tmpdist;
   if ((fLastDistanceToOutWithV.p == p) && (fLastDistanceToOutWithV.vec == v) )
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
   else
   {
     if (currentside == kSurface)
     {
       // particle is just on a boundary.
       // If the particle is exiting from the volume, return 0.
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
         return fLastDistanceToOutWithV.value;
       }
     }
   }
      
   // now, we can take smallest positive distance.
   
   // Initialize
   //
   G4double distance = kInfinity;
       
   // find intersections and choose nearest one.
   //
   G4VTwistSurface* surfaces[6];
   surfaces[0] = fLatterTwisted;
   surfaces[1] = fFormerTwisted;
   surfaces[2] = fInnerHype;
   surfaces[3] = fOuterHype;
   surfaces[4] = fLowerEndcap;
   surfaces[5] = fUpperEndcap;
   
   G4int besti = -1;
   G4ThreeVector xx;
   G4ThreeVector bestxx;
   for (auto i=0; i<6; ++i)
   {
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

   return fLastDistanceToOutWithV.value;
}


//=====================================================================
//* DistanceToOut (p) ----------------------------------------------

G4double G4TwistedTubs::DistanceToOut( const G4ThreeVector& p ) const
{
   // DistanceToOut(p):
   // Calculate distance to surface of shape from `inside', 
   // allowing for tolerance

   //
   // checking last value
   //
   
   G4ThreeVector* tmpp;
   G4double* tmpdist;
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
         G4double      distance = kInfinity;
   
         // find intersections and choose nearest one.
         G4VTwistSurface* surfaces[6];
         surfaces[0] = fLatterTwisted;
         surfaces[1] = fFormerTwisted;
         surfaces[2] = fInnerHype;
         surfaces[3] = fOuterHype;
         surfaces[4] = fLowerEndcap;
         surfaces[5] = fUpperEndcap;

         G4ThreeVector xx;
         G4ThreeVector bestxx;
         for (auto i=0; i<6; ++i)
         {
            G4double tmpdistance = surfaces[i]->DistanceTo(p, xx);
            if (tmpdistance < distance)
            {
               distance = tmpdistance;
               bestxx = xx;
            }
         }
         *tmpdist = distance;
   
         return fLastDistanceToOut.value;
      }
      default :
      {
         G4Exception("G4TwistedTubs::DistanceToOut(p)", "GeomSolids0003",
                     FatalException, "Unknown point location!");
      }
   } // switch end

   return 0.;
}

//=====================================================================
//* StreamInfo --------------------------------------------------------

std::ostream& G4TwistedTubs::StreamInfo(std::ostream& os) const
{
  //
  // Stream object contents to an output stream
  //
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4TwistedTubs\n"
     << " Parameters: \n"
     << "    -ve end Z              : " << fEndZ[0]/mm << " mm \n"
     << "    +ve end Z              : " << fEndZ[1]/mm << " mm \n"
     << "    inner end radius(-ve z): " << fEndInnerRadius[0]/mm << " mm \n"
     << "    inner end radius(+ve z): " << fEndInnerRadius[1]/mm << " mm \n"
     << "    outer end radius(-ve z): " << fEndOuterRadius[0]/mm << " mm \n"
     << "    outer end radius(+ve z): " << fEndOuterRadius[1]/mm << " mm \n"
     << "    inner radius (z=0)     : " << fInnerRadius/mm << " mm \n"
     << "    outer radius (z=0)     : " << fOuterRadius/mm << " mm \n"
     << "    twisted angle          : " << fPhiTwist/degree << " degrees \n"
     << "    inner stereo angle     : " << fInnerStereo/degree << " degrees \n"
     << "    outer stereo angle     : " << fOuterStereo/degree << " degrees \n"
     << "    phi-width of a piece   : " << fDPhi/degree << " degrees \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


//=====================================================================
//* DiscribeYourselfTo ------------------------------------------------

void G4TwistedTubs::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddSolid (*this);
}

//=====================================================================
//* GetExtent ---------------------------------------------------------

G4VisExtent G4TwistedTubs::GetExtent() const
{
  // Define the sides of the box into which the G4Tubs instance would fit.
  //
  G4ThreeVector pmin,pmax;
  BoundingLimits(pmin,pmax);
  return G4VisExtent(pmin.x(),pmax.x(),
                     pmin.y(),pmax.y(),
                     pmin.z(),pmax.z());
}

//=====================================================================
//* CreatePolyhedron --------------------------------------------------

G4Polyhedron* G4TwistedTubs::CreatePolyhedron () const 
{
  // number of meshes
  //
  G4double absPhiTwist = std::abs(fPhiTwist);
  G4double dA = std::max(fDPhi,absPhiTwist);
  const G4int k =
    G4int(G4Polyhedron::GetNumberOfRotationSteps() * dA / twopi) + 2;
  const G4int n =
    G4int(G4Polyhedron::GetNumberOfRotationSteps() * absPhiTwist / twopi) + 2;

  const G4int nnodes = 4*(k-1)*(n-2) + 2*k*k ;
  const G4int nfaces = 4*(k-1)*(n-1) + 2*(k-1)*(k-1) ;

  G4Polyhedron* ph = new G4Polyhedron;
  typedef G4double G4double3[3];
  typedef G4int G4int4[4];
  G4double3* xyz = new G4double3[nnodes];  // number of nodes 
  G4int4*  faces = new G4int4[nfaces] ;    // number of faces
  fLowerEndcap->GetFacets(k,k,xyz,faces,0) ;
  fUpperEndcap->GetFacets(k,k,xyz,faces,1) ;
  fInnerHype->GetFacets(k,n,xyz,faces,2) ;
  fFormerTwisted->GetFacets(k,n,xyz,faces,3) ;
  fOuterHype->GetFacets(k,n,xyz,faces,4) ;
  fLatterTwisted->GetFacets(k,n,xyz,faces,5) ;

  ph->createPolyhedron(nnodes,nfaces,xyz,faces);

  delete[] xyz;
  delete[] faces;

  return ph;
}

//=====================================================================
//* GetPolyhedron -----------------------------------------------------

G4Polyhedron* G4TwistedTubs::GetPolyhedron () const
{
  if (fpPolyhedron == nullptr ||
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
//* CreateSurfaces ----------------------------------------------------

void G4TwistedTubs::CreateSurfaces() 
{
   // create 6 surfaces of TwistedTub

   G4ThreeVector x0(0, 0, fEndZ[0]);
   G4ThreeVector n (0, 0, -1);

   fLowerEndcap = new G4TwistTubsFlatSide("LowerEndcap",
                                    fEndInnerRadius, fEndOuterRadius,
                                    fDPhi, fEndPhi, fEndZ, -1) ;

   fUpperEndcap = new G4TwistTubsFlatSide("UpperEndcap",  
                                    fEndInnerRadius, fEndOuterRadius,
                                    fDPhi, fEndPhi, fEndZ, 1) ;

   G4RotationMatrix    rotHalfDPhi;
   rotHalfDPhi.rotateZ(0.5*fDPhi);

   fLatterTwisted = new G4TwistTubsSide("LatterTwisted",
                                         fEndInnerRadius, fEndOuterRadius,
                                         fDPhi, fEndPhi, fEndZ, 
                                         fInnerRadius, fOuterRadius, fKappa,
                                         1 ) ; 
   fFormerTwisted = new G4TwistTubsSide("FormerTwisted", 
                                         fEndInnerRadius, fEndOuterRadius,
                                         fDPhi, fEndPhi, fEndZ, 
                                         fInnerRadius, fOuterRadius, fKappa,
                                         -1 ) ; 

   fInnerHype = new G4TwistTubsHypeSide("InnerHype",
                                        fEndInnerRadius, fEndOuterRadius,
                                        fDPhi, fEndPhi, fEndZ, 
                                        fInnerRadius, fOuterRadius,fKappa,
                                        fTanInnerStereo, fTanOuterStereo, -1) ;
   fOuterHype = new G4TwistTubsHypeSide("OuterHype", 
                                        fEndInnerRadius, fEndOuterRadius,
                                        fDPhi, fEndPhi, fEndZ, 
                                        fInnerRadius, fOuterRadius,fKappa,
                                        fTanInnerStereo, fTanOuterStereo, 1) ;


   // set neighbour surfaces
   //
   fLowerEndcap->SetNeighbours(fInnerHype, fLatterTwisted,
                               fOuterHype, fFormerTwisted);
   fUpperEndcap->SetNeighbours(fInnerHype, fLatterTwisted,
                               fOuterHype, fFormerTwisted);
   fLatterTwisted->SetNeighbours(fInnerHype, fLowerEndcap,
                                 fOuterHype, fUpperEndcap);
   fFormerTwisted->SetNeighbours(fInnerHype, fLowerEndcap,
                                 fOuterHype, fUpperEndcap);
   fInnerHype->SetNeighbours(fLatterTwisted, fLowerEndcap,
                             fFormerTwisted, fUpperEndcap);
   fOuterHype->SetNeighbours(fLatterTwisted, fLowerEndcap,
                             fFormerTwisted, fUpperEndcap);
}


//=====================================================================
//* GetEntityType -----------------------------------------------------

G4GeometryType G4TwistedTubs::GetEntityType() const
{
  return G4String("G4TwistedTubs");
}

//=====================================================================
//* Clone -------------------------------------------------------------

G4VSolid* G4TwistedTubs::Clone() const
{
  return new G4TwistedTubs(*this);
}

//=====================================================================
//* GetCubicVolume ----------------------------------------------------

G4double G4TwistedTubs::GetCubicVolume()
{
  if (fCubicVolume == 0.)
  {
    G4double DPhi  = GetDPhi();
    G4double Z0    = GetEndZ(0);
    G4double Z1    = GetEndZ(1);
    G4double Ain   = GetInnerRadius();
    G4double Aout  = GetOuterRadius();
    G4double R0in  = GetEndInnerRadius(0);
    G4double R1in  = GetEndInnerRadius(1);
    G4double R0out = GetEndOuterRadius(0);
    G4double R1out = GetEndOuterRadius(1);

    // V_hyperboloid = pi*h*(2*a*a + R*R)/3
    fCubicVolume = (2.*(Z1 - Z0)*(Aout + Ain)*(Aout - Ain)
                    + Z1*(R1out + R1in)*(R1out - R1in)
                    - Z0*(R0out + R0in)*(R0out - R0in))*DPhi/6.;
  }
  return fCubicVolume;
}

//=====================================================================
//* GetLateralArea ----------------------------------------------------

G4double
G4TwistedTubs::GetLateralArea(G4double a, G4double r, G4double z) const
{
  if (z == 0) return 0.;
  G4double h = std::abs(z);
  G4double area = h*a;
  if (std::abs(a - r) > kCarTolerance)
  {
    G4double aa = a*a;
    G4double hh = h*h;
    G4double rr = r*r;
    G4double cc = aa*hh/(rr - aa);
    G4double k  = std::sqrt(aa + cc)/cc;
    G4double kh = k*h;
    area = 0.5*a*(h*std::sqrt(1. + kh*kh) + std::asinh(kh)/k);
  }
  return GetDPhi()*area;
}

//=====================================================================
//* GetPhiCutArea -----------------------------------------------------

G4double
G4TwistedTubs::GetPhiCutArea(G4double a, G4double r, G4double z) const
{
  if (GetDPhi() >= CLHEP::twopi || r <= 0 || z == 0) return 0.;
  G4double h = std::abs(z);
  G4double area = h*a;
  if (GetPhiTwist() > kCarTolerance)
  {
    G4double sinw = std::sin(0.5*GetPhiTwist())*h/GetZHalfLength();
    G4double p = sinw*r/h;
    G4double q = sinw*r/a;
    G4double pp = p*p;
    G4double qq = q*q;
    G4double pq = p*q;
    G4double sqroot = std::sqrt(pp + qq + 1);
    area = (pq*sqroot +
            0.5*p*(pp + 3.)*std::atanh(q/sqroot) +
            0.5*q*(qq + 3.)*std::atanh(p/sqroot) +
            std::atan(sqroot/(pq)) - CLHEP::halfpi)*h*a/(3.*pq);
  }
  return area;
}

//=====================================================================
//* GetSurfaceArea ----------------------------------------------------

G4double G4TwistedTubs::GetSurfaceArea()
{
  if (fSurfaceArea == 0.)
  {
    G4double dphi = GetDPhi();
    G4double Ainn = GetInnerRadius();
    G4double Aout = GetOuterRadius();
    G4double Rinn0 = GetEndInnerRadius(0);
    G4double Rout0 = GetEndOuterRadius(0);
    G4double Rinn1 = GetEndInnerRadius(1);
    G4double Rout1 = GetEndOuterRadius(1);
    G4double z0 = GetEndZ(0);
    G4double z1 = GetEndZ(1);

    G4double base0 = 0.5*dphi*(Rout0*Rout0 - Rinn0*Rinn0); // lower base
    G4double inner0 = GetLateralArea(Ainn, Rinn0, z0); // lower inner surface
    G4double outer0 = GetLateralArea(Aout, Rout0, z0); // lower outer surface
    G4double cut0 =                                    // lower phi cut
      GetPhiCutArea(Aout, Rout0, z0) - GetPhiCutArea(Ainn, Rinn0, z0);

    G4double base1 = base0;
    G4double inner1 = inner0;
    G4double outer1 = outer0;
    G4double cut1 = cut0;
    if (std::abs(z0) != std::abs(z1))
    {
      base1 = 0.5*dphi*(Rout1*Rout1 - Rinn1*Rinn1); // upper base
      inner1 = GetLateralArea(Ainn, Rinn1, z1); // upper inner surface
      outer1 = GetLateralArea(Aout, Rout1, z1); // upper outer surface
      cut1 =                                    // upper phi cut
      GetPhiCutArea(Aout, Rout1, z1) - GetPhiCutArea(Ainn, Rinn1, z1);
    }
    fSurfaceArea = base0 + base1 +
      ((z0*z1 < 0) ?
      (inner0 + inner1 + outer0 + outer1 + 2.*(cut0 + cut1)) :
      std::abs(inner0 - inner1 + outer0 - outer1 + 2.*(cut0 - cut1)));
  }
  return fSurfaceArea;
}

//=====================================================================
//* GetPointOnSurface -------------------------------------------------

G4ThreeVector G4TwistedTubs::GetPointOnSurface() const
{

  G4double z = G4RandFlat::shoot(fEndZ[0],fEndZ[1]);
  G4double phi , phimin, phimax ;
  G4double x   , xmin,   xmax ;
  G4double r   , rmin,   rmax ;

  G4double a1 = fOuterHype->GetSurfaceArea() ;
  G4double a2 = fInnerHype->GetSurfaceArea() ;
  G4double a3 = fLatterTwisted->GetSurfaceArea() ;
  G4double a4 = fFormerTwisted->GetSurfaceArea() ;
  G4double a5 = fLowerEndcap->GetSurfaceArea()  ;
  G4double a6 = fUpperEndcap->GetSurfaceArea() ;

  G4double chose = G4RandFlat::shoot(0.,a1 + a2 + a3 + a4 + a5 + a6) ;

  if(chose < a1)
  {

    phimin = fOuterHype->GetBoundaryMin(z) ;
    phimax = fOuterHype->GetBoundaryMax(z) ;
    phi = G4RandFlat::shoot(phimin,phimax) ;

    return fOuterHype->SurfacePoint(phi,z,true) ;

  }
  else if ( (chose >= a1) && (chose < a1 + a2 ) )
  {

    phimin = fInnerHype->GetBoundaryMin(z) ;
    phimax = fInnerHype->GetBoundaryMax(z) ;
    phi = G4RandFlat::shoot(phimin,phimax) ;

    return fInnerHype->SurfacePoint(phi,z,true) ;

  }
  else if ( (chose >= a1 + a2 ) && (chose < a1 + a2 + a3 ) ) 
  {

    xmin = fLatterTwisted->GetBoundaryMin(z) ; 
    xmax = fLatterTwisted->GetBoundaryMax(z) ;
    x = G4RandFlat::shoot(xmin,xmax) ;
    
    return fLatterTwisted->SurfacePoint(x,z,true) ;

  }
  else if ( (chose >= a1 + a2 + a3  ) && (chose < a1 + a2 + a3 + a4  ) )
  {

    xmin = fFormerTwisted->GetBoundaryMin(z) ; 
    xmax = fFormerTwisted->GetBoundaryMax(z) ;
    x = G4RandFlat::shoot(xmin,xmax) ;

    return fFormerTwisted->SurfacePoint(x,z,true) ;
   }
  else if( (chose >= a1 + a2 + a3 + a4  )&&(chose < a1 + a2 + a3 + a4 + a5 ) )
  {
    rmin = GetEndInnerRadius(0) ;
    rmax = GetEndOuterRadius(0) ;
    r = std::sqrt(G4RandFlat::shoot()*(sqr(rmax)-sqr(rmin))+sqr(rmin));

    phimin = fLowerEndcap->GetBoundaryMin(r) ; 
    phimax = fLowerEndcap->GetBoundaryMax(r) ;
    phi    = G4RandFlat::shoot(phimin,phimax) ;

    return fLowerEndcap->SurfacePoint(phi,r,true) ;
  }
  else
  {
    rmin = GetEndInnerRadius(1) ;
    rmax = GetEndOuterRadius(1) ;
    r = rmin + (rmax-rmin)*std::sqrt(G4RandFlat::shoot());

    phimin = fUpperEndcap->GetBoundaryMin(r) ; 
    phimax = fUpperEndcap->GetBoundaryMax(r) ;
    phi    = G4RandFlat::shoot(phimin,phimax) ;

    return fUpperEndcap->SurfacePoint(phi,r,true) ;
  }
}
