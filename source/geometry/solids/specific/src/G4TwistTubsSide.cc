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
// $Id: G4TwistTubsSide.cc 72937 2013-08-14 13:20:38Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistTubsSide.cc
//
// Author: 
//   01-Aug-2002 - Kotoyo Hoshina (hoshina@hepburn.s.chiba-u.ac.jp)
//
// History:
//   13-Nov-2003 - O.Link (Oliver.Link@cern.ch), Integration in Geant4
//                 from original version in Jupiter-2.5.02 application.
//   29-Apr-2004 - O.Link. Bug fixed in GetAreaCode
// --------------------------------------------------------------------

#include "G4TwistTubsSide.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistTubsSide::G4TwistTubsSide(const G4String         &name,
                                   const G4RotationMatrix &rot,
                                   const G4ThreeVector    &tlate,
                                         G4int             handedness,
                                   const G4double          kappa,
                                   const EAxis             axis0,
                                   const EAxis             axis1,
                                         G4double          axis0min,
                                         G4double          axis1min,
                                         G4double          axis0max,
                                         G4double          axis1max)
   : G4VTwistSurface(name, rot, tlate, handedness, axis0, axis1,
                axis0min, axis1min, axis0max, axis1max),
     fKappa(kappa)
{
   if (axis0 == kZAxis && axis1 == kXAxis) {
      G4Exception("G4TwistTubsSide::G4TwistTubsSide()", "GeomSolids0002",
                  FatalErrorInArgument, "Should swap axis0 and axis1!");
   }
   fIsValidNorm = false;
   SetCorners();
   SetBoundaries();
}

G4TwistTubsSide::G4TwistTubsSide(const G4String     &name,
                                         G4double      EndInnerRadius[2],
                                         G4double      EndOuterRadius[2],
                                         G4double      DPhi,
                                         G4double      EndPhi[2],
                                         G4double      EndZ[2], 
                                         G4double      InnerRadius,
                                         G4double      OuterRadius,
                                         G4double      Kappa,
                                         G4int         handedness)
  : G4VTwistSurface(name)
{  
   fHandedness = handedness;   // +z = +ve, -z = -ve
   fAxis[0]    = kXAxis; // in local coordinate system
   fAxis[1]    = kZAxis;
   fAxisMin[0] = InnerRadius;  // Inner-hype radius at z=0
   fAxisMax[0] = OuterRadius;  // Outer-hype radius at z=0
   fAxisMin[1] = EndZ[0];
   fAxisMax[1] = EndZ[1];

   fKappa = Kappa;
   fRot.rotateZ( fHandedness > 0
                 ? -0.5*DPhi
                 :  0.5*DPhi );
   fTrans.set(0, 0, 0);
   fIsValidNorm = false;
   
   SetCorners( EndInnerRadius, EndOuterRadius, EndPhi, EndZ) ;
   SetBoundaries();
}


//=====================================================================
//* Fake default constructor ------------------------------------------

G4TwistTubsSide::G4TwistTubsSide( __void__& a )
  : G4VTwistSurface(a), fKappa(0.)
{
}


//=====================================================================
//* destructor --------------------------------------------------------

G4TwistTubsSide::~G4TwistTubsSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistTubsSide::GetNormal(const G4ThreeVector &tmpxx, 
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
   
   G4ThreeVector er(1, fKappa * xx.z(), 0);
   G4ThreeVector ez(0, fKappa * xx.x(), 1);
   G4ThreeVector normal = fHandedness*(er.cross(ez));

   if (isGlobal) {
      fCurrentNormal.normal = ComputeGlobalDirection(normal.unit());
   } else {
      fCurrentNormal.normal = normal.unit();
   }
   return fCurrentNormal.normal;
}

//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int G4TwistTubsSide::DistanceToSurface(const G4ThreeVector &gp,
                                          const G4ThreeVector &gv,
                                                G4ThreeVector  gxx[],
                                                G4double       distance[],
                                                G4int          areacode[],
                                                G4bool         isvalid[],
                                                EValidate      validate)
{
   // Coordinate system:
   //
   //    The coordinate system is so chosen that the intersection of
   //    the twisted surface with the z=0 plane coincides with the
   //    x-axis. 
   //    Rotation matrix from this coordinate system (local system)
   //    to global system is saved in fRot field.
   //    So the (global) particle position and (global) velocity vectors, 
   //    p and v, should be rotated fRot.inverse() in order to convert
   //    to local vectors.
   //
   // Equation of a twisted surface:
   //
   //    x(rho(z=0), z) = rho(z=0)
   //    y(rho(z=0), z) = rho(z=0)*K*z
   //    z(rho(z=0), z) = z
   //    with
   //       K = std::tan(fPhiTwist/2)/fZHalfLen
   //
   // Equation of a line:
   //
   //    gxx = p + t*v
   //    with
   //       p = fRot.inverse()*gp  
   //       v = fRot.inverse()*gv
   //
   // Solution for intersection:
   //
   //    Required time for crossing is given by solving the
   //    following quadratic equation:
   //
   //       a*t^2 + b*t + c = 0
   //
   //    where
   //
   //       a = K*v_x*v_z
   //       b = K*(v_x*p_z + v_z*p_x) - v_y
   //       c = K*p_x*p_z - p_y
   //
   //    Out of the possible two solutions you must choose
   //    the one that gives a positive rho(z=0).
   //
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
   // special case! 
   // p is origin or
   //

   G4double absvz = std::fabs(v.z());

   if ((absvz < DBL_MIN) && (std::fabs(p.x() * v.y() - p.y() * v.x()) < DBL_MIN)) {
      // no intersection

      isvalid[0] = false;
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid[0], 0, validate, &gp, &gv);
      return 0;
   } 
   
   // 
   // special case end
   //

   
   G4double a = fKappa * v.x() * v.z();
   G4double b = fKappa * (v.x() * p.z() + v.z() * p.x()) - v.y();
   G4double c = fKappa * p.x() * p.z() - p.y();
   G4double D = b * b - 4 * a * c;             // discriminant
   G4int vout = 0;

   if (std::fabs(a) < DBL_MIN) {
      if (std::fabs(b) > DBL_MIN) { 

         // single solution

         distance[0] = - c / b;
         xx[0]       = p + distance[0]*v;
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
            // we must omit x(rho,z) = rho(z=0) < 0
            if (xx[0].x() > 0) {
               areacode[0] = sInside;
               if (distance[0] >= 0) isvalid[0] = true;
            } else {
               distance[0] = kInfinity;
               fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0],
                                              areacode[0], isvalid[0],
                                              0, validate, &gp, &gv);
               return vout;
            } 
         }
                 
         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 1, validate, &gp, &gv);
         vout = 1;

      } else {
         // if a=b=0 , v.y=0 and (v.x=0 && p.x=0) or (v.z=0 && p.z=0) .
         //    if v.x=0 && p.x=0, no intersection unless p is on z-axis
         //    (in that case, v is paralell to surface). 
         //    if v.z=0 && p.z=0, no intersection unless p is on x-axis
         //    (in that case, v is paralell to surface). 
         // return distance = infinity.

         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 0, validate, &gp, &gv);
      }
      
   } else if (D > DBL_MIN) {   

      // double solutions

      D = std::sqrt(D);
      G4double      factor = 0.5/a;
      G4double      tmpdist[2] = {kInfinity, kInfinity};
      G4ThreeVector tmpxx[2];
      G4int         tmpareacode[2] = {sOutside, sOutside};
      G4bool        tmpisvalid[2]  = {false, false};
      G4int i;

      for (i=0; i<2; i++) {
         G4double bminusD = - b - D;

         // protection against round off error  
         //G4double protection = 1.0e-6;
         G4double protection = 0;
         if ( b * D < 0 && std::fabs(bminusD / D) < protection ) {
            G4double acovbb = (a*c)/(b*b);
            tmpdist[i] = - c/b * ( 1 - acovbb * (1 + 2*acovbb));
         } else { 
            tmpdist[i] = factor * bminusD;
         }

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
            // we must choose x(rho,z) = rho(z=0) > 0
            if (tmpxx[i].x() > 0) {
               tmpareacode[i] = sInside;
               if (tmpdist[i] >= 0) tmpisvalid[i] = true;
               continue;
            } else {
               tmpdist[i] = kInfinity;
               continue;
            }                     
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

      // protection against roundoff error

      for (G4int k=0; k<2; k++) {
         if (!isvalid[k]) continue;

         G4ThreeVector xxonsurface(xx[k].x(), fKappa * std::fabs(xx[k].x())
                                              * xx[k].z() , xx[k].z());
         G4double      deltaY  =  (xx[k] - xxonsurface).mag();

         if ( deltaY > 0.5*kCarTolerance ) {

           G4int maxcount = 10;
           G4int l;
           G4double      lastdeltaY = deltaY; 
           for (l=0; l<maxcount; l++) {
             G4ThreeVector surfacenormal = GetNormal(xxonsurface); 
             distance[k] = DistanceToPlaneWithV(p, v, xxonsurface,
                                                surfacenormal, xx[k]);
             deltaY      = (xx[k] - xxonsurface).mag();
             if (deltaY > lastdeltaY) {
               
             }
             gxx[k]      = ComputeGlobalPoint(xx[k]);

               if (deltaY <= 0.5*kCarTolerance) {
            
                 break;
               }
               xxonsurface.set(xx[k].x(),
                               fKappa * std::fabs(xx[k].x()) * xx[k].z(),
                               xx[k].z());
            }
            if (l == maxcount) {
               std::ostringstream message;
               message << "Exceeded maxloop count!" << G4endl
                      << "        maxloop count " << maxcount;
               G4Exception("G4TwistTubsFlatSide::DistanceToSurface(p,v)",
                           "GeomSolids0003",  FatalException, message);
            }
         }

      }
      vout = 2;
   } else {
      // if D<0, no solution
      // if D=0, just grazing the surfaces, return kInfinity

      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                     isvalid[0], 0, validate, &gp, &gv);
   }

   return vout;
}

//=====================================================================
//* DistanceToSurface -------------------------------------------------

G4int G4TwistTubsSide::DistanceToSurface(const G4ThreeVector &gp,
                                                G4ThreeVector  gxx[],
                                                G4double       distance[],
                                                G4int          areacode[])
{  
   fCurStat.ResetfDone(kDontValidate, &gp);
   G4int i = 0;
   if (fCurStat.IsDone()) {
      for (i=0; i<fCurStat.GetNXX(); i++) {
         gxx[i] = fCurStat.GetXX(i);
         distance[i] = fCurStat.GetDistance(i);
         areacode[i] = fCurStat.GetAreacode(i);
      }
      return fCurStat.GetNXX();
   } else {
      // initialize
      for (i=0; i<2; i++) {
         distance[i] = kInfinity;
         areacode[i] = sOutside;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   
   const G4double halftol = 0.5 * kCarTolerance; 

   G4ThreeVector  p       = ComputeLocalPoint(gp);
   G4ThreeVector  xx;
   G4int          parity  = (fKappa >= 0 ? 1 : -1);
 
   // 
   // special case! 
   // If p is on surface, or
   // p is on z-axis, 
   // return here immediatery.
   //
   
   G4ThreeVector  lastgxx[2];
   for (i=0; i<2; i++) {
      lastgxx[i] = fCurStatWithV.GetXX(i);
   } 

   if  ((gp - lastgxx[0]).mag() < halftol
     || (gp - lastgxx[1]).mag() < halftol) { 
      // last winner, or last poststep point is on the surface.
      xx = p;
      distance[0] = 0;
      gxx[0] = gp;

      G4bool isvalid = true;
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                             isvalid, 1, kDontValidate, &gp);
      return 1;
   }
          
   if (p.getRho() == 0) { 
      // p is on z-axis. Namely, p is on twisted surface (invalid area).
      // We must return here, however, returning distance to x-minimum
      // boundary is better than return 0-distance.
      //
      G4bool isvalid = true;
      if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
         distance[0] = DistanceToBoundary(sAxis0 & sAxisMin, xx, p);
         areacode[0] = sInside;
      } else {
         distance[0] = 0;
         xx.set(0., 0., 0.);
      }
      gxx[0] = ComputeGlobalPoint(xx);
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 0, kDontValidate, &gp);
      return 1;
   } 

   // 
   // special case end
   //

   // set corner points of quadrangle try area ...

   G4ThreeVector A;  // foot of normal from p to boundary of sAxis0 & sAxisMin
   G4ThreeVector C;  // foot of normal from p to boundary of sAxis0 & sAxisMax
   G4ThreeVector B;       // point on boundary sAxis0 & sAxisMax at z = A.z()
   G4ThreeVector D;       // point on boundary sAxis0 & sAxisMin at z = C.z()

   // G4double      distToA; // distance from p to A
   DistanceToBoundary(sAxis0 & sAxisMin, A, p);
   // G4double      distToC; // distance from p to C 
   DistanceToBoundary(sAxis0 & sAxisMax, C, p);
   
   // is p.z between a.z and c.z?
   // p.z must be bracketed a.z and c.z.
   if (A.z() > C.z()) {
      if (p.z() > A.z()) {
         A = GetBoundaryAtPZ(sAxis0 & sAxisMin, p);
      } else if (p.z() < C.z()) {
         C = GetBoundaryAtPZ(sAxis0 & sAxisMax, p);
      }
   } else {
      if (p.z() > C.z()) {
         C = GetBoundaryAtPZ(sAxis0 & sAxisMax, p);
      } else if (p.z() < A.z()) {
         A = GetBoundaryAtPZ(sAxis0 & sAxisMin, p);
      }
   }

   G4ThreeVector  d[2];     // direction vectors of boundary
   G4ThreeVector  x0[2];    // foot of normal from line to p 
   G4int          btype[2]; // boundary type

   for (i=0; i<2; i++) {
      if (i == 0) {
         GetBoundaryParameters((sAxis0 & sAxisMax), d[i], x0[i], btype[i]);
         B = x0[i] + ((A.z() - x0[i].z()) / d[i].z()) * d[i]; 
         // x0 + t*d , d is direction unit vector.
      } else {
         GetBoundaryParameters((sAxis0 & sAxisMin), d[i], x0[i], btype[i]);
         D = x0[i] + ((C.z() - x0[i].z()) / d[i].z()) * d[i]; 
      }
   }

   // In order to set correct diagonal, swap A and D, C and B if needed.  
   G4ThreeVector pt(p.x(), p.y(), 0.);
   G4double      rc = std::fabs(p.x());
   G4ThreeVector surfacevector(rc, rc * fKappa * p.z(), 0.); 
   G4int         pside = AmIOnLeftSide(pt, surfacevector); 
   G4double      test  = (A.z() - C.z()) * parity * pside;  

   if (test == 0) {
      if (pside == 0) {
         // p is on surface.
         xx = p;
         distance[0] = 0;
         gxx[0] = gp;

         G4bool isvalid = true;
         fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 1, kDontValidate, &gp);
         return 1;
      } else {
         // A.z = C.z(). return distance to line.
         d[0] = C - A;
         distance[0] = DistanceToLine(p, A, d[0], xx);
         areacode[0] = sInside;
         gxx[0] = ComputeGlobalPoint(xx);
         G4bool isvalid = true;
         fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 1, kDontValidate, &gp);
         return 1;
      } 

   } else if (test < 0) {

      // wrong diagonal. vector AC is crossing the surface!  
      // swap A and D, C and B
      G4ThreeVector tmp;
      tmp = A;
      A = D;
      D = tmp;
      tmp = C;
      C = B;
      B = tmp; 

   } else {
      // correct diagonal. nothing to do.  
   }

   // Now, we chose correct diaglnal.
   // First try. divide quadrangle into double triangle by diagonal and 
   // calculate distance to both surfaces.

   G4ThreeVector xxacb;   // foot of normal from plane ACB to p
   G4ThreeVector nacb;    // normal of plane ACD
   G4ThreeVector xxcad;   // foot of normal from plane CAD to p
   G4ThreeVector ncad;    // normal of plane CAD
   G4ThreeVector AB(A.x(), A.y(), 0);
   G4ThreeVector DC(C.x(), C.y(), 0);

   G4double distToACB = G4VTwistSurface::DistanceToPlane(p, A, C-A, AB, xxacb, nacb) * parity;
   G4double distToCAD = G4VTwistSurface::DistanceToPlane(p, C, C-A, DC, xxcad, ncad) * parity;

   // if calculated distance = 0, return  

   if (std::fabs(distToACB) <= halftol || std::fabs(distToCAD) <= halftol) {
      xx = (std::fabs(distToACB) < std::fabs(distToCAD) ? xxacb : xxcad); 
      areacode[0] = sInside;
      gxx[0] = ComputeGlobalPoint(xx);
      distance[0] = 0;
      G4bool isvalid = true;
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0] , areacode[0],
                                isvalid, 1, kDontValidate, &gp);
      return 1;
   }
   
   if (distToACB * distToCAD > 0 && distToACB < 0) {
      // both distToACB and distToCAD are negative.
      // divide quadrangle into double triangle by diagonal
      G4ThreeVector normal;
      distance[0] = DistanceToPlane(p, A, B, C, D, parity, xx, normal);
   } else {
      if (distToACB * distToCAD > 0) {
         // both distToACB and distToCAD are positive.
         // Take smaller one.
         if (distToACB <= distToCAD) {
            distance[0] = distToACB;
            xx   = xxacb;
         } else {
            distance[0] = distToCAD;
            xx   = xxcad;
         }
      } else {
         // distToACB * distToCAD is negative.
         // take positive one
         if (distToACB > 0) {
            distance[0] = distToACB;
            xx   = xxacb;
         } else {
            distance[0] = distToCAD;
            xx   = xxcad;
         }
      }
      
   }
   areacode[0] = sInside;
   gxx[0]      = ComputeGlobalPoint(xx);
   G4bool isvalid = true;
   fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                             isvalid, 1, kDontValidate, &gp);
   return 1;
}

//=====================================================================
//* DistanceToPlane ---------------------------------------------------

G4double G4TwistTubsSide::DistanceToPlane(const G4ThreeVector &p,
                                           const G4ThreeVector &A,
                                           const G4ThreeVector &B,
                                           const G4ThreeVector &C,
                                           const G4ThreeVector &D,
                                           const G4int          parity,
                                                 G4ThreeVector &xx,
                                                 G4ThreeVector &n)
{
   const G4double halftol = 0.5 * kCarTolerance;
   
   G4ThreeVector M = 0.5*(A + B);
   G4ThreeVector N = 0.5*(C + D);
   G4ThreeVector xxanm;  // foot of normal from p to plane ANM
   G4ThreeVector nanm;   // normal of plane ANM
   G4ThreeVector xxcmn;  // foot of normal from p to plane CMN
   G4ThreeVector ncmn;   // normal of plane CMN

   G4double distToanm = G4VTwistSurface::DistanceToPlane(p, A, (N - A), (M - A), xxanm, nanm) * parity;
   G4double distTocmn = G4VTwistSurface::DistanceToPlane(p, C, (M - C), (N - C), xxcmn, ncmn) * parity;

   // if p is behind of both surfaces, abort.
   if (distToanm * distTocmn > 0 && distToanm < 0) {
     G4Exception("G4TwistTubsSide::DistanceToPlane()",
                 "GeomSolids0003", FatalException,
                 "Point p is behind the surfaces.");
   }

   // if p is on surface, return 0.
   if (std::fabs(distToanm) <= halftol) {
      xx = xxanm;
      n  = nanm * parity;
      return 0;
   } else if (std::fabs(distTocmn) <= halftol) {
      xx = xxcmn;
      n  = ncmn * parity;
      return 0;
   }
   
   if (distToanm <= distTocmn) {
      if (distToanm > 0) {
         // both distanses are positive. take smaller one.
         xx = xxanm;
         n  = nanm * parity;
         return distToanm;
      } else {
         // take -ve distance and call the function recursively.
         return DistanceToPlane(p, A, M, N, D, parity, xx, n);
      }
   } else {
      if (distTocmn > 0) {
         // both distanses are positive. take smaller one.
         xx = xxcmn;
         n  = ncmn * parity;
         return distTocmn;
      } else {
         // take -ve distance and call the function recursively.
         return DistanceToPlane(p, C, N, M, B, parity, xx, n);
      }
   }
}

//=====================================================================
//* GetAreaCode -------------------------------------------------------

G4int G4TwistTubsSide::GetAreaCode(const G4ThreeVector &xx, 
                                          G4bool withTol)
{
   // We must use the function in local coordinate system.
   // See the description of DistanceToSurface(p,v).
   
   const G4double ctol = 0.5 * kCarTolerance;
   G4int areacode = sInside;
   
   if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
      G4int xaxis = 0;
      G4int zaxis = 1;
      
      if (withTol) {

         G4bool isoutside   = false;

         // test boundary of xaxis

         if (xx.x() < fAxisMin[xaxis] + ctol) {
            areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary; 
            if (xx.x() <= fAxisMin[xaxis] - ctol) isoutside = true;

         } else if (xx.x() > fAxisMax[xaxis] - ctol) {
            areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
            if (xx.x() >= fAxisMax[xaxis] + ctol)  isoutside = true;
         }

         // test boundary of z-axis

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

         // if isoutside = true, clear inside bit.             
         // if not on boundary, add axis information.             
         
         if (isoutside) {
            G4int tmpareacode = areacode & (~sInside);
            areacode = tmpareacode;
         } else if ((areacode & sBoundary) != sBoundary) {
            areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisZ);
         }           
         
      } else {

         // boundary of x-axis

         if (xx.x() < fAxisMin[xaxis] ) {
            areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary;
         } else if (xx.x() > fAxisMax[xaxis]) {
            areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
         }
         
         // boundary of z-axis

         if (xx.z() < fAxisMin[zaxis]) {
            areacode |= (sAxis1 & (sAxisZ | sAxisMin));
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary; 
           
         } else if (xx.z() > fAxisMax[zaxis]) {
            areacode |= (sAxis1 & (sAxisZ | sAxisMax)) ;
            if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
            else                        areacode |= sBoundary; 
         }

         if ((areacode & sBoundary) != sBoundary) {
            areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisZ);
         }           
      }
      return areacode;
   } else {
      G4Exception("G4TwistTubsSide::GetAreaCode()",
                  "GeomSolids0001", FatalException,
                  "Feature NOT implemented !");
   }
   return areacode;
}

//=====================================================================
//* SetCorners( arglist ) -------------------------------------------------

void G4TwistTubsSide::SetCorners(
                                  G4double      endInnerRad[2],
                                  G4double      endOuterRad[2],
                                  G4double      endPhi[2],
                                  G4double      endZ[2])
{
   // Set Corner points in local coodinate.   

   if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
   
      G4int zmin = 0 ;  // at -ve z
      G4int zmax = 1 ;  // at +ve z

      G4double x, y, z;
      
      // corner of Axis0min and Axis1min
      x = endInnerRad[zmin]*std::cos(endPhi[zmin]);
      y = endInnerRad[zmin]*std::sin(endPhi[zmin]);
      z = endZ[zmin];
      SetCorner(sC0Min1Min, x, y, z);
      
      // corner of Axis0max and Axis1min
      x = endOuterRad[zmin]*std::cos(endPhi[zmin]);
      y = endOuterRad[zmin]*std::sin(endPhi[zmin]);
      z = endZ[zmin];
      SetCorner(sC0Max1Min, x, y, z);
      
      // corner of Axis0max and Axis1max
      x = endOuterRad[zmax]*std::cos(endPhi[zmax]);
      y = endOuterRad[zmax]*std::sin(endPhi[zmax]);
      z = endZ[zmax];
      SetCorner(sC0Max1Max, x, y, z);
      
      // corner of Axis0min and Axis1max
      x = endInnerRad[zmax]*std::cos(endPhi[zmax]);
      y = endInnerRad[zmax]*std::sin(endPhi[zmax]);
      z = endZ[zmax];
      SetCorner(sC0Min1Max, x, y, z);

   } else {
      std::ostringstream message;
      message << "Feature NOT implemented !" << G4endl
              << "        fAxis[0] = " << fAxis[0] << G4endl
              << "        fAxis[1] = " << fAxis[1];
      G4Exception("G4TwistTubsSide::SetCorners()",
                  "GeomSolids0001", FatalException, message);
   }
}

//=====================================================================
//* SetCorners() ------------------------------------------------------

void G4TwistTubsSide::SetCorners()
{
   G4Exception("G4TwistTubsSide::SetCorners()",
               "GeomSolids0001", FatalException,
               "Method NOT implemented !");
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistTubsSide::SetBoundaries()
{
   // Set direction-unit vector of boundary-lines in local coodinate. 
   //   
   G4ThreeVector direction;
   
   if (fAxis[0] == kXAxis && fAxis[1] == kZAxis) {
      
      // sAxis0 & sAxisMin
      direction = GetCorner(sC0Min1Max) - GetCorner(sC0Min1Min);
      direction = direction.unit();
      SetBoundary(sAxis0 & (sAxisX | sAxisMin), direction, 
                  GetCorner(sC0Min1Min), sAxisZ) ;
      
      // sAxis0 & sAxisMax
      direction = GetCorner(sC0Max1Max) - GetCorner(sC0Max1Min);
      direction = direction.unit();
      SetBoundary(sAxis0 & (sAxisX | sAxisMax), direction, 
                  GetCorner(sC0Max1Min), sAxisZ);
                  
      // sAxis1 & sAxisMin
      direction = GetCorner(sC0Max1Min) - GetCorner(sC0Min1Min);
      direction = direction.unit();
      SetBoundary(sAxis1 & (sAxisZ | sAxisMin), direction, 
                  GetCorner(sC0Min1Min), sAxisX);
                  
      // sAxis1 & sAxisMax
      direction = GetCorner(sC0Max1Max) - GetCorner(sC0Min1Max);
      direction = direction.unit();
      SetBoundary(sAxis1 & (sAxisZ | sAxisMax), direction, 
                  GetCorner(sC0Min1Max), sAxisX);
                  
   } else {
      std::ostringstream message;
      message << "Feature NOT implemented !" << G4endl
              << "        fAxis[0] = " << fAxis[0] << G4endl
              << "        fAxis[1] = " << fAxis[1];
      G4Exception("G4TwistTubsSide::SetCorners()",
                  "GeomSolids0001", FatalException, message);
   }
}

//=====================================================================
//* GetFacets() -------------------------------------------------------

void G4TwistTubsSide::GetFacets( G4int k, G4int n, G4double xyz[][3],
                                 G4int faces[][4], G4int iside ) 
{

  G4double z ;     // the two parameters for the surface equation
  G4double x,xmin,xmax ;

  G4ThreeVector p ;  // a point on the surface, given by (z,u)

  G4int nnode ;
  G4int nface ;

  // calculate the (n-1)*(k-1) vertices

  G4int i,j ;

  for ( i = 0 ; i<n ; i++ )
  {

    z = fAxisMin[1] + i*(fAxisMax[1]-fAxisMin[1])/(n-1) ;

    for ( j = 0 ; j<k ; j++ ) {

      nnode = GetNode(i,j,k,n,iside) ;

      xmin = GetBoundaryMin(z) ; 
      xmax = GetBoundaryMax(z) ;

      if (fHandedness < 0) { 
        x = xmin + j*(xmax-xmin)/(k-1) ;
      } else {               
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
