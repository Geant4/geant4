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
// $Id: G4TwistTubsHypeSide.cc 72937 2013-08-14 13:20:38Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistTubsHypeSide.cc
//
// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
// --------------------------------------------------------------------

#include "G4TwistTubsHypeSide.hh"
#include "G4PhysicalConstants.hh"
#include "G4GeometryTolerance.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistTubsHypeSide::G4TwistTubsHypeSide(const G4String         &name,
                                         const G4RotationMatrix &rot,
                                         const G4ThreeVector    &tlate,
                                         const G4int             handedness,
                                         const G4double          kappa,
                                         const G4double          tanstereo,
                                         const G4double          r0,
                                         const EAxis             axis0,
                                         const EAxis             axis1,
                                               G4double          axis0min,
                                               G4double          axis1min,
                                               G4double          axis0max,
                                               G4double          axis1max )
  : G4VTwistSurface(name, rot, tlate, handedness, axis0, axis1,
                   axis0min, axis1min, axis0max, axis1max),
    fKappa(kappa), fTanStereo(tanstereo),
    fTan2Stereo(tanstereo*tanstereo), fR0(r0), fR02(r0*r0), fDPhi(twopi)
{
   if ( (axis0 == kZAxis) && (axis1 == kPhi) )
   {
      G4Exception("G4TwistTubsHypeSide::G4TwistTubsHypeSide()",
                  "GeomSolids0002", FatalErrorInArgument,
                  "Should swap axis0 and axis1!");
   }
   
   fInside.gp.set(kInfinity, kInfinity, kInfinity);
   fInside.inside = kOutside;
   fIsValidNorm = false;
   
   SetCorners();
   SetBoundaries();

}

G4TwistTubsHypeSide::G4TwistTubsHypeSide(const G4String      &name,
                                         G4double         EndInnerRadius[2],
                                         G4double         EndOuterRadius[2],
                                         G4double         DPhi,
                                         G4double         EndPhi[2],
                                         G4double         EndZ[2], 
                                         G4double         InnerRadius,
                                         G4double         OuterRadius,
                                         G4double         Kappa,
                                         G4double         TanInnerStereo,
                                         G4double         TanOuterStereo,
                                         G4int            handedness)
   : G4VTwistSurface(name)
{

   fHandedness = handedness;   // +z = +ve, -z = -ve
   fAxis[0]    = kPhi;
   fAxis[1]    = kZAxis;
   fAxisMin[0] = kInfinity;         // we cannot fix boundary min of Phi, 
   fAxisMax[0] = kInfinity;         // because it depends on z.
   fAxisMin[1] = EndZ[0];
   fAxisMax[1] = EndZ[1];
   fKappa      = Kappa;
   fDPhi       = DPhi ;

   if (handedness < 0) { // inner hyperbolic surface
      fTanStereo  = TanInnerStereo;
      fR0         = InnerRadius;
   } else {              // outer hyperbolic surface
      fTanStereo  = TanOuterStereo;
      fR0         = OuterRadius;
   }
   fTan2Stereo = fTanStereo * fTanStereo;
   fR02        = fR0 * fR0;
   
   fTrans.set(0, 0, 0);
   fIsValidNorm = false;

   fInside.gp.set(kInfinity, kInfinity, kInfinity);
   fInside.inside = kOutside;
   
   SetCorners(EndInnerRadius, EndOuterRadius, DPhi, EndPhi, EndZ) ; 

   SetBoundaries();
}

//=====================================================================
//* Fake default constructor ------------------------------------------

G4TwistTubsHypeSide::G4TwistTubsHypeSide( __void__& a )
  : G4VTwistSurface(a), fKappa(0.), fTanStereo(0.), fTan2Stereo(0.),
    fR0(0.), fR02(0.), fDPhi(0.)
{
}

//=====================================================================
//* destructor --------------------------------------------------------

G4TwistTubsHypeSide::~G4TwistTubsHypeSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistTubsHypeSide::GetNormal(const G4ThreeVector &tmpxx, 
                                                   G4bool isGlobal) 
{
   // GetNormal returns a normal vector at a surface (or very close
   // to surface) point at tmpxx.
   // If isGlobal=true, it returns the normal in global coordinate.
   //
   
   G4ThreeVector xx;
   if (isGlobal) {
      xx = ComputeLocalPoint(tmpxx);
      if ((xx - fCurrentNormal.p).mag() < 0.5 * kCarTolerance) {
         return ComputeGlobalDirection(fCurrentNormal.normal);
      }
   } else {
      xx = tmpxx;
      if (xx == fCurrentNormal.p) {
         return fCurrentNormal.normal;
      }
   }
   
   fCurrentNormal.p = xx;

   G4ThreeVector normal( xx.x(), xx.y(), -xx.z() * fTan2Stereo);
   normal *= fHandedness;
   normal = normal.unit();

   if (isGlobal) {
      fCurrentNormal.normal = ComputeLocalDirection(normal);
   } else {
      fCurrentNormal.normal = normal;
   }
   return fCurrentNormal.normal;
}

//=====================================================================
//* Inside() ----------------------------------------------------------

EInside G4TwistTubsHypeSide::Inside(const G4ThreeVector &gp) 
{
   // Inside returns 
   const G4double halftol
     = 0.5 * G4GeometryTolerance::GetInstance()->GetRadialTolerance();

   if (fInside.gp == gp) {
      return fInside.inside;
   }
   fInside.gp = gp;
   
   G4ThreeVector p = ComputeLocalPoint(gp);
  

   if (p.mag() < DBL_MIN) {
      fInside.inside = kOutside;
      return fInside.inside;
   }
   
   G4double rhohype = GetRhoAtPZ(p);
   G4double distanceToOut = fHandedness * (rhohype - p.getRho());
                            // +ve : inside

   if (distanceToOut < -halftol) {

     fInside.inside = kOutside;

   } else {

      G4int areacode = GetAreaCode(p);
      if (IsOutside(areacode)) {
         fInside.inside = kOutside;
      } else if (IsBoundary(areacode)) {
         fInside.inside = kSurface;
      } else if (IsInside(areacode)) {
         if (distanceToOut <= halftol) {
            fInside.inside = kSurface;
         } else {
            fInside.inside = kInside;
         }
      } else {
         G4cout << "WARNING - G4TwistTubsHypeSide::Inside()" << G4endl
                << "          Invalid option !" << G4endl
                << "          name, areacode, distanceToOut = "
                << GetName() << ", " << std::hex << areacode << std::dec << ", "
                << distanceToOut << G4endl;
      }
   }
   
   return fInside.inside; 
}

//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int G4TwistTubsHypeSide::DistanceToSurface(const G4ThreeVector &gp,
                                             const G4ThreeVector &gv,
                                                   G4ThreeVector  gxx[],
                                                   G4double       distance[],
                                                   G4int          areacode[],
                                                   G4bool         isvalid[],
                                                   EValidate      validate)
{
   //
   // Decide if and where a line intersects with a hyperbolic
   // surface (of infinite extent)
   //
   // Arguments:
   //     p       - (in) Point on trajectory
   //     v       - (in) Vector along trajectory
   //     r2      - (in) Square of radius at z = 0
   //     tan2phi - (in) std::tan(stereo)**2
   //     s       - (out) Up to two points of intersection, where the
   //                     intersection point is p + s*v, and if there are
   //                     two intersections, s[0] < s[1]. May be negative.
   // Returns:
   //     The number of intersections. If 0, the trajectory misses.
   //
   //
   // Equation of a line:
   //
   //       x = x0 + s*tx      y = y0 + s*ty      z = z0 + s*tz
   //
   // Equation of a hyperbolic surface:
   //
   //       x**2 + y**2 = r**2 + (z*tanPhi)**2
   //
   // Solution is quadratic:
   //
   //  a*s**2 + b*s + c = 0
   //
   // where:
   //
   //  a = tx**2 + ty**2 - (tz*tanPhi)**2
   //
   //  b = 2*( x0*tx + y0*ty - z0*tz*tanPhi**2 )
   //
   //  c = x0**2 + y0**2 - r**2 - (z0*tanPhi)**2
   //
      
   fCurStatWithV.ResetfDone(validate, &gp, &gv);

   if (fCurStatWithV.IsDone()) {
      G4int i;
      for (i=0; i<fCurStatWithV.GetNXX(); i++) {
         gxx[i] = fCurStatWithV.GetXX(i);
         distance[i] = fCurStatWithV.GetDistance(i);
         areacode[i] = fCurStatWithV.GetAreacode(i);
         isvalid[i]  = fCurStatWithV.IsValid(i);
      }
      return fCurStatWithV.GetNXX();
   } else {
      // initialize
      G4int i;
      for (i=0; i<2; i++) {
         distance[i] = kInfinity;
         areacode[i] = sOutside;
         isvalid[i]  = false;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   
   G4ThreeVector p = ComputeLocalPoint(gp);
   G4ThreeVector v = ComputeLocalDirection(gv);
   G4ThreeVector xx[2]; 
   
   //
   // special case!  p is on origin.
   // 

   if (p.mag() == 0) {
      // p is origin. 
      // unique solution of 2-dimension question in r-z plane 
      // Equations:
      //    r^2 = fR02 + z^2*fTan2Stere0
      //    r = beta*z
      //        where 
      //        beta = vrho / vz
      // Solution (z value of intersection point):
      //    xxz = +- std::sqrt (fR02 / (beta^2 - fTan2Stereo))
      //

      G4double vz    = v.z();
      G4double absvz = std::fabs(vz);
      G4double vrho  = v.getRho();       
      G4double vslope = vrho/vz;
      G4double vslope2 = vslope * vslope;
      if (vrho == 0 || (vrho/absvz) <= (absvz*std::fabs(fTanStereo)/absvz)) {
         // vz/vrho is bigger than slope of asymptonic line
         distance[0] = kInfinity;
         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 0, validate, &gp, &gv);
         return 0;
      }
       
      if (vz) { 
         G4double xxz  = std::sqrt(fR02 / (vslope2 - fTan2Stereo)) 
                        * (vz / std::fabs(vz)) ;
         G4double t = xxz / vz;
         xx[0].set(t*v.x(), t*v.y(), xxz);
      } else {
         // p.z = 0 && v.z =0
         xx[0].set(v.x()*fR0, v.y()*fR0, 0);  // v is a unit vector.
      }
      distance[0] = xx[0].mag();
      gxx[0]      = ComputeGlobalPoint(xx[0]);

      if (validate == kValidateWithTol) {
         areacode[0] = GetAreaCode(xx[0]);
         if (!IsOutside(areacode[0])) {
            if (distance[0] >= 0) isvalid[0] = true;
         }
      } else if (validate == kValidateWithoutTol) {
         areacode[0] = GetAreaCode(xx[0], false);
         if (IsInside(areacode[0])) {
            if (distance[0] >= 0) isvalid[0] = true;
         }
      } else { // kDontValidate                       
         areacode[0] = sInside;
            if (distance[0] >= 0) isvalid[0] = true;
      }
                 
      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 1, validate, &gp, &gv);
      return 1;
   }

   //
   // special case end.
   // 

   G4double a = v.x()*v.x() + v.y()*v.y() - v.z()*v.z()*fTan2Stereo;
   G4double b = 2.0 * ( p.x() * v.x() + p.y() * v.y() - p.z() * v.z() * fTan2Stereo );
   G4double c = p.x()*p.x() + p.y()*p.y() - fR02 - p.z()*p.z()*fTan2Stereo;
   G4double D = b*b - 4*a*c;          //discriminant
   G4int vout = 0;
   
   if (std::fabs(a) < DBL_MIN) {
      if (std::fabs(b) > DBL_MIN) {           // single solution

         distance[0] = -c/b;
         xx[0] = p + distance[0]*v;
         gxx[0] = ComputeGlobalPoint(xx[0]);

         if (validate == kValidateWithTol) {
            areacode[0] = GetAreaCode(xx[0]);
            if (!IsOutside(areacode[0])) {
               if (distance[0] >= 0) isvalid[0] = true;
            }
         } else if (validate == kValidateWithoutTol) {
            areacode[0] = GetAreaCode(xx[0], false);
            if (IsInside(areacode[0])) {
               if (distance[0] >= 0) isvalid[0] = true;
            }
         } else { // kDontValidate                       
            areacode[0] = sInside;
               if (distance[0] >= 0) isvalid[0] = true;
         }
                 
         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 1, validate, &gp, &gv);
         vout = 1;
         
      } else {
         // if a=b=0 and c != 0, p is origin and v is parallel to asymptotic line.
         // if a=b=c=0, p is on surface and v is paralell to stereo wire. 
         // return distance = infinity.

         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 0, validate, &gp, &gv);

         vout = 0;
      }
      
   } else if (D > DBL_MIN) {         // double solutions
      
      D = std::sqrt(D);
      G4double      factor = 0.5/a;
      G4double      tmpdist[2] = {kInfinity, kInfinity};
      G4ThreeVector tmpxx[2] ;
      G4int         tmpareacode[2] = {sOutside, sOutside};
      G4bool        tmpisvalid[2]  = {false, false};
      G4int i;

      for (i=0; i<2; i++) {
         tmpdist[i] = factor*(-b - D);
         D = -D;
         tmpxx[i] = p + tmpdist[i]*v;
        
         if (validate == kValidateWithTol) {
            tmpareacode[i] = GetAreaCode(tmpxx[i]);
            if (!IsOutside(tmpareacode[i])) {
               if (tmpdist[i] >= 0) tmpisvalid[i] = true;
               continue;
            }
         } else if (validate == kValidateWithoutTol) {
            tmpareacode[i] = GetAreaCode(tmpxx[i], false);
            if (IsInside(tmpareacode[i])) {
               if (tmpdist[i] >= 0) tmpisvalid[i] = true;
               continue;
            }
         } else { // kDontValidate
            tmpareacode[i] = sInside;
               if (tmpdist[i] >= 0) tmpisvalid[i] = true;
            continue;
         }
      }      

      if (tmpdist[0] <= tmpdist[1]) {
          distance[0] = tmpdist[0];
          distance[1] = tmpdist[1];
          xx[0]       = tmpxx[0];
          xx[1]       = tmpxx[1];
          gxx[0]      = ComputeGlobalPoint(tmpxx[0]);
          gxx[1]      = ComputeGlobalPoint(tmpxx[1]);
          areacode[0] = tmpareacode[0];
          areacode[1] = tmpareacode[1];
          isvalid[0]  = tmpisvalid[0];
          isvalid[1]  = tmpisvalid[1];
      } else {
          distance[0] = tmpdist[1];
          distance[1] = tmpdist[0];
          xx[0]       = tmpxx[1];
          xx[1]       = tmpxx[0];
          gxx[0]      = ComputeGlobalPoint(tmpxx[1]);
          gxx[1]      = ComputeGlobalPoint(tmpxx[0]);
          areacode[0] = tmpareacode[1];
          areacode[1] = tmpareacode[0];
          isvalid[0]  = tmpisvalid[1];
          isvalid[1]  = tmpisvalid[0];
      }
         
      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                     isvalid[0], 2, validate, &gp, &gv);
      fCurStatWithV.SetCurrentStatus(1, gxx[1], distance[1], areacode[1],
                                     isvalid[1], 2, validate, &gp, &gv);
      vout = 2;
      
   } else {
      // if D<0, no solution
      // if D=0, just grazing the surfaces, return kInfinity

      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                     isvalid[0], 0, validate, &gp, &gv);
      vout = 0;
   }
   return vout;
}

   
//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int G4TwistTubsHypeSide::DistanceToSurface(const G4ThreeVector &gp,
                                                   G4ThreeVector  gxx[],
                                                   G4double       distance[],
                                                   G4int          areacode[])
{
    // Find the approximate distance of a point of a hyperbolic surface.
    // The distance must be an underestimate. 
    // It will also be nice (although not necessary) that the estimate is
    // always finite no matter how close the point is.
    //
    // We arranged G4Hype::ApproxDistOutside and G4Hype::ApproxDistInside
    // for this function. See these discriptions.
    
   const G4double halftol
     = 0.5 * G4GeometryTolerance::GetInstance()->GetRadialTolerance();

   fCurStat.ResetfDone(kDontValidate, &gp);

   if (fCurStat.IsDone()) {
      for (G4int i=0; i<fCurStat.GetNXX(); i++) {
         gxx[i] = fCurStat.GetXX(i);
         distance[i] = fCurStat.GetDistance(i);
         areacode[i] = fCurStat.GetAreacode(i);
      }
      return fCurStat.GetNXX();
   } else {
      // initialize
      for (G4int i=0; i<2; i++) {
         distance[i] = kInfinity;
         areacode[i] = sOutside;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   

   G4ThreeVector p = ComputeLocalPoint(gp);
   G4ThreeVector xx;

   //
   // special case!
   // If p is on surface, return distance = 0 immediatery .
   //
   G4ThreeVector  lastgxx[2];
   for (G4int i=0; i<2; i++) {
      lastgxx[i] = fCurStatWithV.GetXX(i);
   }

   if ((gp - lastgxx[0]).mag() < halftol || (gp - lastgxx[1]).mag() < halftol) {
      // last winner, or last poststep point is on the surface.
      xx = p;             
      gxx[0] = gp;
      distance[0] = 0;      

      G4bool isvalid = true;
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 1, kDontValidate, &gp);

      return 1;

   }
   //
   // special case end
   //
       
   G4double prho       = p.getRho();
   G4double pz         = std::fabs(p.z());           // use symmetry
   G4double r1         = std::sqrt(fR02 + pz * pz * fTan2Stereo);
   
   G4ThreeVector pabsz(p.x(), p.y(), pz);
    
   if (prho > r1 + halftol) {  // p is outside of Hyperbolic surface

      // First point xx1
      G4double t = r1 / prho;
      G4ThreeVector xx1(t * pabsz.x(), t * pabsz.y() , pz);
      
      // Second point xx2
      G4double z2 = (prho * fTanStereo + pz) / (1 + fTan2Stereo);
      G4double r2 = std::sqrt(fR02 + z2 * z2 * fTan2Stereo);
      t = r2 / prho;
      G4ThreeVector xx2(t * pabsz.x(), t * pabsz.y() , z2);
            
      G4double len = (xx2 - xx1).mag();
      if (len < DBL_MIN) {
         // xx2 = xx1?? I guess we
         // must have really bracketed the normal
         distance[0] = (pabsz - xx1).mag();
         xx = xx1;
      } else {
         distance[0] = DistanceToLine(pabsz, xx1, (xx2 - xx1) , xx);
      }
      
   } else if (prho < r1 - halftol) { // p is inside of Hyperbolic surface.
           
      // First point xx1
      G4double t;
      G4ThreeVector xx1;
      if (prho < DBL_MIN) {
         xx1.set(r1, 0. , pz);
      } else {
         t = r1 / prho;
         xx1.set(t * pabsz.x(), t * pabsz.y() , pz);
      }
      
      // dr, dz is tangential vector of Hyparbolic surface at xx1
      // dr = r, dz = z*tan2stereo
      G4double dr  = pz * fTan2Stereo;
      G4double dz  = r1;
      G4double tanbeta   = dr / dz;
      G4double pztanbeta = pz * tanbeta;
      
      // Second point xx2 
      // xx2 is intersection between x-axis and tangential vector
      G4double r2 = r1 - pztanbeta;
      G4ThreeVector xx2;
      if (prho < DBL_MIN) {
         xx2.set(r2, 0. , 0.);
      } else {
         t  = r2 / prho;
         xx2.set(t * pabsz.x(), t * pabsz.y() , 0.);
      }
      
      G4ThreeVector d = xx2 - xx1;
      distance[0] = DistanceToLine(pabsz, xx1, d, xx);
          
   } else {  // p is on Hyperbolic surface.
   
      distance[0] = 0;
      xx.set(p.x(), p.y(), pz);

   }

   if (p.z() < 0) {
      G4ThreeVector tmpxx(xx.x(), xx.y(), -xx.z());
      xx = tmpxx;
   }
       
   gxx[0] = ComputeGlobalPoint(xx);
   areacode[0]    = sInside;
   G4bool isvalid = true;
   fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                             isvalid, 1, kDontValidate, &gp);
   return 1;
}

//=====================================================================
//* GetAreaCode -------------------------------------------------------

G4int G4TwistTubsHypeSide::GetAreaCode(const G4ThreeVector &xx, 
                                             G4bool         withTol)
{
   const G4double ctol = 0.5 * kCarTolerance;
   G4int areacode = sInside;

   if ((fAxis[0] == kPhi && fAxis[1] == kZAxis))  {
      //G4int phiaxis = 0;
      G4int zaxis   = 1;
      
      if (withTol) {

         G4bool isoutside      = false;
         G4int  phiareacode    = GetAreaCodeInPhi(xx);
         G4bool isoutsideinphi = IsOutside(phiareacode);

         // test boundary of phiaxis

         if ((phiareacode & sAxisMin) == sAxisMin) {

            areacode |= (sAxis0 & (sAxisPhi | sAxisMin)) | sBoundary;
            if (isoutsideinphi) isoutside = true;

         } else if ((phiareacode & sAxisMax)  == sAxisMax) {

            areacode |= (sAxis0 & (sAxisPhi | sAxisMax)) | sBoundary;
            if (isoutsideinphi) isoutside = true;

         }

         // test boundary of zaxis

         if (xx.z() < fAxisMin[zaxis] + ctol) {

            areacode |= (sAxis1 & (sAxisZ | sAxisMin));
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary;

            if (xx.z() <= fAxisMin[zaxis] - ctol) isoutside = true;

         } else if (xx.z() > fAxisMax[zaxis] - ctol) {

            areacode |= (sAxis1 & (sAxisZ | sAxisMax));
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary;

            if (xx.z() >= fAxisMax[zaxis] + ctol) isoutside = true;
         }

         // if isoutside = true, clear sInside bit.
         // if not on boundary, add boundary information. 

         if (isoutside) {
            G4int tmpareacode = areacode & (~sInside);
            areacode = tmpareacode;
         } else if ((areacode & sBoundary) != sBoundary) {
            areacode |= (sAxis0 & sAxisPhi) | (sAxis1 & sAxisZ);
         }

         return areacode;
      
      } else {

         G4int phiareacode = GetAreaCodeInPhi(xx, false);
         
         // test boundary of z-axis

         if (xx.z() < fAxisMin[zaxis]) {

            areacode |= (sAxis1 & (sAxisZ | sAxisMin)) | sBoundary;

         } else if (xx.z() > fAxisMax[zaxis]) {

            areacode |= (sAxis1 & (sAxisZ | sAxisMax)) | sBoundary;

         }

         // boundary of phi-axis

         if (phiareacode == sAxisMin) {

            areacode |= (sAxis0 & (sAxisPhi | sAxisMin));
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary; 
             
         } else if (phiareacode == sAxisMax) {

            areacode |= (sAxis0 & (sAxisPhi | sAxisMax));
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary; 
           
         }

         // if not on boundary, add boundary information. 

         if ((areacode & sBoundary) != sBoundary) {
            areacode |= (sAxis0 & sAxisPhi) | (sAxis1 & sAxisZ);
         }
         return areacode;
      }
   } else {
      std::ostringstream message;
      message << "Feature NOT implemented !" << G4endl
              << "        fAxis[0] = " << fAxis[0] << G4endl
              << "        fAxis[1] = " << fAxis[1];
      G4Exception("G4TwistTubsHypeSide::GetAreaCode()",
                  "GeomSolids0001", FatalException, message);
   }
   return areacode;
}

//=====================================================================
//* GetAreaCodeInPhi --------------------------------------------------

G4int G4TwistTubsHypeSide::GetAreaCodeInPhi(const G4ThreeVector &xx,
                                                  G4bool withTol)
{
   
   G4ThreeVector lowerlimit; // lower phi-boundary limit at z = xx.z()
   G4ThreeVector upperlimit; // upper phi-boundary limit at z = xx.z()
   lowerlimit = GetBoundaryAtPZ(sAxis0 & sAxisMin, xx);
   upperlimit = GetBoundaryAtPZ(sAxis0 & sAxisMax, xx);

   G4int  areacode  = sInside;
   G4bool isoutside = false; 
   
   if (withTol) {
         
      if (AmIOnLeftSide(xx, lowerlimit) >= 0) {        // xx is on lowerlimit
         areacode |= (sAxisMin | sBoundary);
         if (AmIOnLeftSide(xx, lowerlimit) > 0) isoutside = true; 

      } else if (AmIOnLeftSide(xx, upperlimit) <= 0) { // xx is on upperlimit
         areacode |= (sAxisMax | sBoundary);
         if (AmIOnLeftSide(xx, upperlimit) < 0) isoutside = true; 
      }

      // if isoutside = true, clear inside bit.

      if (isoutside) {
         G4int tmpareacode = areacode & (~sInside);
         areacode = tmpareacode;
      }


   } else {
   
      if (AmIOnLeftSide(xx, lowerlimit, false) >= 0) {
         areacode |= (sAxisMin | sBoundary);
      } else if (AmIOnLeftSide(xx, upperlimit, false) <= 0) {
         areacode |= (sAxisMax | sBoundary);
      }
   }

   return areacode;
   
}

//=====================================================================
//* SetCorners(EndInnerRadius, EndOuterRadius,DPhi,EndPhi,EndZ) -------

void G4TwistTubsHypeSide::SetCorners(
                                     G4double         EndInnerRadius[2],
                                     G4double         EndOuterRadius[2],
                                     G4double         DPhi,
                                     G4double         endPhi[2],
                                     G4double         endZ[2] 
                                     )
{
   // Set Corner points in local coodinate.

   if (fAxis[0] == kPhi && fAxis[1] == kZAxis) {

      G4int i;
      G4double endRad[2];
      G4double halfdphi = 0.5*DPhi;
      
      for (i=0; i<2; i++) { // i=0,1 : -ve z, +ve z
         endRad[i] = (fHandedness == 1 ? EndOuterRadius[i]
                                      : EndInnerRadius[i]);
      }

      G4int zmin = 0 ;  // at -ve z
      G4int zmax = 1 ;  // at +ve z

      G4double x, y, z;
      
      // corner of Axis0min and Axis1min
      x = endRad[zmin]*std::cos(endPhi[zmin] - halfdphi);
      y = endRad[zmin]*std::sin(endPhi[zmin] - halfdphi);
      z = endZ[zmin];
      SetCorner(sC0Min1Min, x, y, z);
      
      // corner of Axis0max and Axis1min
      x = endRad[zmin]*std::cos(endPhi[zmin] + halfdphi);
      y = endRad[zmin]*std::sin(endPhi[zmin] + halfdphi);
      z = endZ[zmin];
      SetCorner(sC0Max1Min, x, y, z);
      
      // corner of Axis0max and Axis1max
      x = endRad[zmax]*std::cos(endPhi[zmax] + halfdphi);
      y = endRad[zmax]*std::sin(endPhi[zmax] + halfdphi);
      z = endZ[zmax];
      SetCorner(sC0Max1Max, x, y, z);
      
      // corner of Axis0min and Axis1max
      x = endRad[zmax]*std::cos(endPhi[zmax] - halfdphi);
      y = endRad[zmax]*std::sin(endPhi[zmax] - halfdphi);
      z = endZ[zmax];
      SetCorner(sC0Min1Max, x, y, z);

   } else {
      std::ostringstream message;
      message << "Feature NOT implemented !" << G4endl
              << "        fAxis[0] = " << fAxis[0] << G4endl
              << "        fAxis[1] = " << fAxis[1];
      G4Exception("G4TwistTubsHypeSide::SetCorners()",
                  "GeomSolids0001", FatalException, message);
   }
}


//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistTubsHypeSide::SetCorners()
{
   G4Exception("G4TwistTubsHypeSide::SetCorners()",
               "GeomSolids0001", FatalException,
               "Method NOT implemented !");
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistTubsHypeSide::SetBoundaries()
{
   // Set direction-unit vector of phi-boundary-lines in local coodinate.
   // sAxis0 must be kPhi.
   // This fanction set lower phi-boundary and upper phi-boundary.

   if (fAxis[0] == kPhi && fAxis[1] == kZAxis) {

      G4ThreeVector direction;
      // sAxis0 & sAxisMin
      direction = GetCorner(sC0Min1Max) - GetCorner(sC0Min1Min);
      direction = direction.unit();
      SetBoundary(sAxis0 & (sAxisPhi | sAxisMin), direction, 
                   GetCorner(sC0Min1Min), sAxisZ);

      // sAxis0 & sAxisMax
      direction = GetCorner(sC0Max1Max) - GetCorner(sC0Max1Min);
      direction = direction.unit();
      SetBoundary(sAxis0 & (sAxisPhi | sAxisMax), direction, 
                  GetCorner(sC0Max1Min), sAxisZ);

      // sAxis1 & sAxisMin
      direction = GetCorner(sC0Max1Min) - GetCorner(sC0Min1Min);
      direction = direction.unit();
      SetBoundary(sAxis1 & (sAxisZ | sAxisMin), direction, 
                   GetCorner(sC0Min1Min), sAxisPhi);

      // sAxis1 & sAxisMax
      direction = GetCorner(sC0Max1Max) - GetCorner(sC0Min1Max);
      direction = direction.unit();
      SetBoundary(sAxis1 & (sAxisZ | sAxisMax), direction, 
                  GetCorner(sC0Min1Max), sAxisPhi);
   } else {
      std::ostringstream message;
      message << "Feature NOT implemented !" << G4endl
              << "        fAxis[0] = " << fAxis[0] << G4endl
              << "        fAxis[1] = " << fAxis[1];
      G4Exception("G4TwistTubsHypeSide::SetBoundaries()",
                  "GeomSolids0001", FatalException, message);
   }
}

//=====================================================================
//* GetFacets() -------------------------------------------------------

void G4TwistTubsHypeSide::GetFacets( G4int k, G4int n, G4double xyz[][3],
                                     G4int faces[][4], G4int iside ) 
{

  G4double z ;     // the two parameters for the surface equation
  G4double x,xmin,xmax ;

  G4ThreeVector p ;  // a point on the surface, given by (z,u)

  G4int nnode ;
  G4int nface ;

  // calculate the (n-1)*(k-1) vertices

  G4int i,j ;

  for ( i = 0 ; i<n ; i++ ) {

    z = fAxisMin[1] + i*(fAxisMax[1]-fAxisMin[1])/(n-1) ;

    for ( j = 0 ; j<k ; j++ )
    {
      nnode = GetNode(i,j,k,n,iside) ;

      xmin = GetBoundaryMin(z) ; 
      xmax = GetBoundaryMax(z) ;

      if (fHandedness < 0) { // inner hyperbolic surface
        x = xmin + j*(xmax-xmin)/(k-1) ;
      } else {               // outer hyperbolic surface
        x = xmax - j*(xmax-xmin)/(k-1) ;
      }

      p = SurfacePoint(x,z,true) ;  // surface point in global coord.system

      xyz[nnode][0] = p.x() ;
      xyz[nnode][1] = p.y() ;
      xyz[nnode][2] = p.z() ;

      if ( i<n-1 && j<k-1 ) {   // clock wise filling
        
        nface = GetFace(i,j,k,n,iside) ;

	faces[nface][0] = GetEdgeVisibility(i,j,k,n,0,1) * ( GetNode(i  ,j  ,k,n,iside)+1) ;  
	faces[nface][1] = GetEdgeVisibility(i,j,k,n,1,1) * ( GetNode(i+1,j  ,k,n,iside)+1) ;
	faces[nface][2] = GetEdgeVisibility(i,j,k,n,2,1) * ( GetNode(i+1,j+1,k,n,iside)+1) ;
	faces[nface][3] = GetEdgeVisibility(i,j,k,n,3,1) * ( GetNode(i  ,j+1,k,n,iside)+1) ;

      }
    }
  }
}
