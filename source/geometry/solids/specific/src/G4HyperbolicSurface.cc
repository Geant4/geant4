/*
 *  G4HyperbolicSurface.cc
 *  
 *
 *  Created by Kotoyo Hoshina on Thu Aug 01 2002.
 *  Copyright (c) 2001 __MyCompanyName__. All rights reserved.
 *
 */

//#define __SOLIDDEBUG__
//#define __SOLIDDEBUGAREACODE__

#include "G4HyperbolicSurface.hh"
#include "G4TwistedTubs.hh"

//=====================================================================
//* constructor -------------------------------------------------------

G4HyperbolicSurface::G4HyperbolicSurface(const G4String         &name,
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
                    :G4VSurface(name, rot, tlate, handedness, axis0, axis1,
                                axis0min, axis1min, axis0max, axis1max),
                     fKappa(kappa), fTanStereo(tanstereo),
                     fTan2Stereo(tanstereo*tanstereo), fR0(r0), fR02(r0*r0)
{
   if (axis0 == kZAxis && axis1 == kPhi) {
      G4cerr << "G4HyperbolicSurface::Constructor: " 
             << "Swap axis0 and axis1. abort. " << G4endl;
      abort();
   }
   
   fInside.gp.set(kInfinity, kInfinity, kInfinity);
   fInside.inside = ::kOutside;
   fIsValidNorm = FALSE;
   
   SetCorners();
   SetBoundaries();

}

G4HyperbolicSurface::G4HyperbolicSurface(const G4String      &name,
                                               G4TwistedTubs *solid,
                                               G4int          handedness)
                    :G4VSurface(name, solid)
{

   fHandedness = handedness;   // +z = +ve, -z = -ve
   fAxis[0]    = kPhi;
   fAxis[1]    = kZAxis;
   fAxisMin[0] = kInfinity;         // we cannot fix boundary min of Phi, 
   fAxisMax[0] = kInfinity;         // because it depends on z.
   fAxisMin[1] = solid->GetEndZ(0);
   fAxisMax[1] = solid->GetEndZ(1);
   fKappa      = solid->GetKappa();

   if (handedness < 0) { // inner hyperbolic surface
      fTanStereo  = solid->GetTanInnerStereo();
      fR0         = solid->GetInnerRadius();
   } else {              // outer hyperbolic surface
      fTanStereo  = solid->GetTanOuterStereo();
      fR0         = solid->GetOuterRadius();
   }
   fTan2Stereo = fTanStereo * fTanStereo;
   fR02        = fR0 * fR0;
   
   fTrans.set(0, 0, 0);
   fIsValidNorm = FALSE;

   fInside.gp.set(kInfinity, kInfinity, kInfinity);
   fInside.inside = ::kOutside;
   
   SetCorners(solid);
   SetBoundaries();
}

//=====================================================================
//* destructor --------------------------------------------------------
G4HyperbolicSurface::~G4HyperbolicSurface()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------
G4ThreeVector G4HyperbolicSurface::GetNormal(const G4ThreeVector &tmpxx, 
                                                   G4bool isGlobal) 
{
   // GetNormal returns a normal vector at a surface (or very close
   // to surface) point at tmpxx.
   // If isGlobal=TRUE, it returns the normal in global coordinate.
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
EInside G4HyperbolicSurface::Inside(const G4ThreeVector &gp) 
{
   // Inside returns 
   static const G4double halftol = 0.5 * kRadTolerance;

   if (fInside.gp == gp) {
      return fInside.inside;
   }
   fInside.gp = gp;
   
   G4ThreeVector p = ComputeLocalPoint(gp);
   
#ifdef __SOLIDDEBUG__
   G4cerr << "      ~~~~~ G4HyperbolicSurface:Inside(gp) start~~~~~~~~~"
          << G4endl;
   G4cerr << "         p : " << p << G4endl;
   G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" 
          << G4endl;
#endif

   if (p.mag() < DBL_MIN) {
      fInside.inside = ::kOutside;
      return fInside.inside;
   }
   
   G4double rhohype = GetRhoAtPZ(p);
   G4double distanceToOut = fHandedness * (rhohype - p.getRho());
                            // +ve : inside

   if (distanceToOut < -halftol) {
#ifdef __SOLIDDEBUG__
      G4cerr << "      G4HyperbolicSurface:Inside(gp) distanceToOut"
             << " is -ve. kOutside. "
             << "distanceToOut = " << distanceToOut << G4endl;
#endif
      fInside.inside = ::kOutside;
   } else {
      G4int areacode = GetAreaCode(p);
      if (IsOutside(areacode)) {
         fInside.inside = ::kOutside;
      } else if (IsBoundary(areacode)) {
         fInside.inside = ::kSurface;
      } else if (IsInside(areacode)) {
         if (distanceToOut <= halftol) {
            fInside.inside = ::kSurface;
         } else {
            fInside.inside = ::kInside;
         }
      } else {
         G4cerr << "      G4HyperbolicSurface::Inside "
                << "invalid option! name, areacode, distanceToOut = "
                << GetName() << " " << hex << areacode << dec << " "
                << distanceToOut << G4endl;
      }
   }
   
#ifdef __SOLIDDEBUG__
   G4cerr << "      ~~~~~ G4HyperbolicSurface:Inside return ~~~~~~~~~~~~~~"
      << G4endl;
   G4cerr << "         Name : " << GetName() << G4endl;
   G4cerr << "         distanceToOut : " << distanceToOut << G4endl;
   G4cerr << "         areacode      : " << std::hex << fInside.inside 
                                         << std::dec << G4endl;
   G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          << G4endl;
#endif
   
   return fInside.inside; 
}

//=====================================================================
//* DistanceToSurface -------------------------------------------------
G4int G4HyperbolicSurface::DistanceToSurface(const G4ThreeVector &gp,
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
   //     tan2phi - (in) tan(stereo)**2
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
         areacode[i] = kOutside;
         isvalid[i]  = FALSE;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   
   G4ThreeVector p = ComputeLocalPoint(gp);
   G4ThreeVector v = ComputeLocalDirection(gv);
   G4ThreeVector xx[2]; 

#ifdef __SOLIDDEBUG__
   G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p,v):Start"
          << G4endl;
   G4cerr << "         Name : " << GetName() << G4endl;
   G4cerr << "         gp   : " << gp << G4endl;
   G4cerr << "         gv   : " << gv << G4endl;
   G4cerr << "         p    : " << p << G4endl;
   G4cerr << "         v    : " << v << G4endl;
   G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          << G4endl;
#endif
   
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
      //    xxz = +- sqrt (fR02 / (beta^2 - fTan2Stereo))
      //

      G4double vz    = v.z();
      G4double absvz = abs(vz);
      G4double vrho  = v.getRho();       
      G4double vslope = vrho/vz;
      G4double vslope2 = vslope * vslope;
      if (vrho == 0 || (vrho/absvz) <= (absvz*fabs(fTanStereo)/absvz)) {
         // vz/vrho is bigger than slope of asymptonic line
         distance[0] = kInfinity;
         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 0, validate, &gp, &gv);
#ifdef __SOLIDDEBUG__
         G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p,v):"
                <<              "return" << G4endl;
         G4cerr << "         vz/vrho is bigger than slope of asymptonic "
                <<           "line." << G4endl;
         G4cerr << "         NAME     : " << GetName() << G4endl;
         G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                <<        "~~~~~" << G4endl;
#endif
         return 0;
      }
       
      if (vz) { 
         G4double xxz  = sqrt(fR02 / (vslope2 - fTan2Stereo)) 
                        * (vz / fabs(vz)) ;
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
            if (distance[0] >= 0) isvalid[0] = TRUE;
         }
      } else if (validate == kValidateWithoutTol) {
         areacode[0] = GetAreaCode(xx[0], FALSE);
         if (IsInside(areacode[0])) {
            if (distance[0] >= 0) isvalid[0] = TRUE;
         }
      } else { // kDontValidate                       
         areacode[0] = kInside;
            if (distance[0] >= 0) isvalid[0] = TRUE;
      }
                 
      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 1, validate, &gp, &gv);
#ifdef __SOLIDDEBUG__
      G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p,v):return"
             << G4endl;
      G4cerr << "         p is on origin. " << G4endl;
      G4cerr << "         NAME        : " << GetName() << G4endl;
      G4cerr << "         xx[0]       : " << xx[0] << G4endl;
      G4cerr << "         gxx[0]      : " << gxx[0] << G4endl;
      G4cerr << "         dist[0]     : " << distance[0] << G4endl;
      G4cerr << "         areacode[0] : " << hex << areacode[0] << G4endl;
      G4cerr << "         isvalid[0]  : " << dec << isvalid[0] << G4endl;
      G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
             << G4endl;
      if (isvalid[0] && GetSolid()->Inside(gxx[0]) != ::kSurface) {
         G4cerr << " valid return value is not on surface! abort." << G4endl;
         abort();
      } 

#endif

      return 1;

   }

   //
   // special case end.
   // 
   

   G4double a = v.x()*v.x() + v.y()*v.y() - v.z()*v.z()*fTan2Stereo;
   G4double b = 2.0 * ( p.x() * v.x() + p.y() * v.y() - p.z() * v.z() * fTan2Stereo );
   G4double c = p.x()*p.x() + p.y()*p.y() - fR02 - p.z()*p.z()*fTan2Stereo;
   G4double D = b*b - 4*a*c;          //discriminant
   
#ifdef __SOLIDDEBUG__
   G4cerr << "      ~~ G4HyperbolicSurface::DistanceToSurface: a,b,c,D ~~~" 
          << G4endl;
   G4cerr << "      //*   NAME      : " << GetName() << G4endl;
   G4cerr << "      //*   p, v      : " << p << " , " << v << G4endl;
   G4cerr << "      //*   a,b,c,D   : " << a << " , " << b << " , " 
                                     << c << " , " << D << G4endl; 
   G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          << G4endl;
#endif  
   
   if (fabs(a) < DBL_MIN) {
      if (fabs(b) > DBL_MIN) {           // single solution

         distance[0] = -c/b;
         xx[0] = p + distance[0]*v;
         gxx[0] = ComputeGlobalPoint(xx[0]);

         if (validate == kValidateWithTol) {
            areacode[0] = GetAreaCode(xx[0]);
            if (!IsOutside(areacode[0])) {
               if (distance[0] >= 0) isvalid[0] = TRUE;
            }
         } else if (validate == kValidateWithoutTol) {
            areacode[0] = GetAreaCode(xx[0], FALSE);
            if (IsInside(areacode[0])) {
               if (distance[0] >= 0) isvalid[0] = TRUE;
            }
         } else { // kDontValidate                       
            areacode[0] = kInside;
               if (distance[0] >= 0) isvalid[0] = TRUE;
         }
                 
         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 1, validate, &gp, &gv);
#ifdef __SOLIDDEBUG__
         G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p,v):return~"
                << G4endl;
         G4cerr << "         Single solution. " << G4endl;
         G4cerr << "         NAME        : " << GetName() << G4endl;
         G4cerr << "         xx[0]       : " << xx[0] << G4endl;
         G4cerr << "         gxx[0]      : " << gxx[0] << G4endl;
         G4cerr << "         dist[0]     : " << distance[0] << G4endl;
         G4cerr << "         areacode[0] : " << hex << areacode[0] << G4endl;
         G4cerr << "         isvalid[0]  : " << dec << isvalid[0] << G4endl;
         G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                << G4endl;
        if (isvalid[0] && GetSolid()->Inside(gxx[0]) != ::kSurface) {
           G4cerr << " valid return value is not on surface! abort." << G4endl;
           abort();
        } 
#endif
         return 1;
         
      } else {
         // if a=b=0 and c != 0, p is origin and v is parallel to asymptotic line.
         // if a=b=c=0, p is on surface and v is paralell to stereo wire. 
         // return distance = infinity.

         fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                        isvalid[0], 0, validate, &gp, &gv);

#ifdef __SOLIDDEBUG__
         G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p,v):return~"
                << G4endl;
         G4cerr << "         No intersection." << G4endl;
         G4cerr << "         a, b, c  : " <<  a  << " , " << b << " , " 
                                          << c << G4endl;
         G4cerr << "         NAME     : " << GetName() << G4endl;
         G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
                << G4endl;
#endif
         return 0;
      }
      
   } else if (D > DBL_MIN) {         // double solutions
      
      D = sqrt(D);
      G4double      factor = 0.5/a;
      G4double      tmpdist[2] = {kInfinity, kInfinity};
      G4ThreeVector tmpxx[2] ;
      G4int         tmpareacode[2] = {kOutside, kOutside};
      G4bool        tmpisvalid[2]  = {FALSE, FALSE};
      G4int i;

      for (i=0; i<2; i++) {
         tmpdist[i] = factor*(-b - D);
         D = -D;
         tmpxx[i] = p + tmpdist[i]*v;
        
         if (validate == kValidateWithTol) {
            tmpareacode[i] = GetAreaCode(tmpxx[i]);
            if (!IsOutside(tmpareacode[i])) {
               if (tmpdist[i] >= 0) tmpisvalid[i] = TRUE;
               continue;
            }
         } else if (validate == kValidateWithoutTol) {
            tmpareacode[i] = GetAreaCode(tmpxx[i], FALSE);
            if (IsInside(tmpareacode[i])) {
               if (tmpdist[i] >= 0) tmpisvalid[i] = TRUE;
               continue;
            }
         } else { // kDontValidate
            tmpareacode[i] = kInside;
               if (tmpdist[i] >= 0) tmpisvalid[i] = TRUE;
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
#ifdef __SOLIDDEBUG__
      G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p,v):return~"
             << G4endl;
      G4cerr << "         NAME,        : " << GetName() << " , " << i << G4endl;
      G4cerr << "         xx[0,1]      : " << xx[0] << " , " << xx[1] << G4endl;
      G4cerr << "         gxx[0,1]     : " << gxx[0] << " , " << gxx[1] << G4endl;
      G4cerr << "         dist[0,1]    : " << distance[0] << " , " << distance[1] << G4endl;
      G4cerr << "         areacode[0,1]: " << hex << areacode[0] << " , " 
                                           << areacode[1] << dec << G4endl;
      G4cerr << "         isvalid[0,1] : " << isvalid[0] << " , " << isvalid[1] 
                                           << G4endl;
      G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
             << G4endl;

      for (G4int k=0; k<2; k++) {
         if (isvalid[k] && GetSolid()->Inside(gxx[k]) != ::kSurface) {
            G4cerr << " valid return value is not on surface! abort. k="
                   << k << G4endl;
            G4ThreeVector pp = gxx[k];
            return DistanceToSurface(pp, gv, gxx, distance, areacode, isvalid, validate); 
            //abort();
         }
      }
#endif
      return 2;
      
   } else {
      // if D<0, no solution
      // if D=0, just grazing the surfaces, return kInfinity

      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                     isvalid[0], 0, validate, &gp, &gv);
#ifdef __SOLIDDEBUG__
      G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p,v):return~"
             << G4endl;
      G4cerr << "         paralell to the surface or on surface but flying"
             <<          "away opposit direction. return 0. " << G4endl;
      G4cerr << "         a, b, c  : " <<  a  << " , " << b << " , " 
                                        << c << G4endl;
      G4cerr << "         NAME     : " << GetName() << G4endl;
      G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
             << G4endl;
#endif
      return 0;
   }
    G4cerr << "      G4HyperbolicSurface::DistanceToSurface(p,v) illigal "
           << "operation!! abort"
           << G4endl;
    abort();          
}

   
//=====================================================================
//* DistanceToSurface -------------------------------------------------
G4int G4HyperbolicSurface::DistanceToSurface(const G4ThreeVector &gp,
                                                   G4ThreeVector  gxx[],
                                                   G4double       distance[],
                                                   G4int          areacode[])
{
    // Find the approximate distance of a point of a hyperbolic surface.
    // The distance must be an underestimate. 
    // It will also be nice (although not necesary) that the estimate is
    // always finite no matter how close the point is.
    //
    // We arranged G4Hype::ApproxDistOutside and G4Hype::ApproxDistInside
    // for this function. See these discriptions.
    
   static const G4double halftol    = 0.5 * kRadTolerance;

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
         areacode[i] = kOutside;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   

   G4ThreeVector p = ComputeLocalPoint(gp);
   G4ThreeVector xx;

#ifdef __SOLIDDEBUG__
   G4cerr << "      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p):Start"
          << G4endl;
   G4cerr << "         Name : " << GetName() << G4endl;
   G4cerr << "         gp   : " << gp << G4endl;
   G4cerr << "         p    : " << p << G4endl;
   G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          << G4endl;
#endif

   //
   // special case!
   // If p is on surface, return distance = 0 immediatery .
   //
   G4ThreeVector  lastgxx[2];
   G4double       distfromlast[2];
   for (G4int i=0; i<2; i++) {
      lastgxx[i] = fCurStatWithV.GetXX(i);
      distfromlast[i] = (gp - lastgxx[i]).mag();
   }

   if ((gp - lastgxx[0]).mag() < halftol || (gp - lastgxx[1]).mag() < halftol) {
      // last winner, or last poststep point is on the surface.
      xx = p;             
      gxx[0] = gp;
      distance[0] = 0;      

#ifdef __BOUNDARYCHECK__
      areacode[0] = GetAreaCode(xx, FALSE);
      if (IsInside(areacode[0])) {
         distance[0] = 0;      
         gxx[0] = gp;
      } else {
         // xx is out of boundary or corner
         if (IsCorner(areacode[0])) {
            xx = GetCorner(areacode[0]);
            distance[0] = (xx - p).mag();
         } else {
            distance[0] = DistanceToBoundary(areacode[0], xx, p);
         }
         gxx[0] = ComputeGlobalPoint(xx);
      }
#endif

      G4bool isvalid = TRUE;
      fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                isvalid, 1, kDontValidate, &gp);
                             
#ifdef __SOLIDDEBUG__
      G4cerr <<"      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p):return"
             << G4endl;
      G4cerr <<"         I'm a last winner ! " << G4endl;
      G4cerr <<"         Otherwise last poststep point is on my surface. "
             << G4endl;
      G4cerr <<"         NAME        : " << GetName() << G4endl;
      G4cerr <<"         xx          : " << xx << G4endl;
      G4cerr <<"         gxx[0]      : " << gxx[0] << G4endl;
      G4cerr <<"         dist[0]     : " << distance[0] << G4endl;
      G4cerr <<"      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
             << G4endl;
#endif

      return 1;

   }
   //
   // special case end
   //
       
   G4double prho       = p.getRho();
   G4double pz         = fabs(p.z());           // use symmetry
   G4double r1         = sqrt(fR02 + pz * pz * fTan2Stereo);
   
   G4ThreeVector pabsz(p.x(), p.y(), pz);
   
#ifdef __SOLIDDEBUG__
   G4cerr << "      ~~~ G4HyperbolicSurface::DistanceToSurface(p):~~~~~~"
          << G4endl;
   G4cerr << "      //*   NAME     : " <<  GetName() << G4endl;
   G4cerr << "      //*   fR02     : " <<  fR02 << G4endl;
   G4cerr << "      //*   p.rho    : " <<  prho << G4endl;
   G4cerr << "      //*   pz       : " <<  pz << G4endl;
   G4cerr << "      //*   r1(z=p.z): " <<  r1 << G4endl;
   G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          << G4endl;
#endif
    
   if (prho > r1 + halftol) {  // p is outside of Hyperbolic surface

      // First point xx1
      G4double t = r1 / prho;
      G4ThreeVector xx1(t * pabsz.x(), t * pabsz.y() , pz);
      
      // Second point xx2
      G4double z2 = (prho * fTanStereo + pz) / (1 + fTan2Stereo);
      G4double r2 = sqrt(fR02 + z2 * z2 * fTan2Stereo);
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
      
#ifdef __SOLIDDEBUG__
      G4cerr << "      ~~~ G4HyperbolicSurface::DistanceToSurface(p):~~~~~~" 
             << G4endl;
      G4cerr << "      //*   p is outside of Hyperbolic surface."  << G4endl;
      G4cerr << "      //*   NAME     : " <<  GetName() << G4endl;
      G4cerr << "      //*   xx1      : " <<  xx1 << G4endl;
      G4cerr << "      //*   xx2      : " <<  xx2 << G4endl;
      G4cerr << "      //*   xx       : " <<  xx  << G4endl;
      G4cerr << "      //*   Len      : " <<  len << G4endl;
      G4cerr << "      //*   Distance : " <<  distance[0] << G4endl;
      G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
             << G4endl;
#endif
            
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
      
#ifdef __SOLIDDEBUG__
      G4cerr << "      ~~~ G4HyperbolicSurface::DistanceToSurface(p):~~~~~~" << G4endl;
      G4cerr << "      //*   p is inside of Hyperbolic surface."  << G4endl;
      G4cerr << "      //*   NAME     : " <<  GetName() << G4endl;
      G4cerr << "      //*   xx1      : " <<  xx1 << G4endl;
      G4cerr << "      //*   xx2      : " <<  xx2 << G4endl;
      G4cerr << "      //*   xx       : " <<  xx  << G4endl;
      G4cerr << "      //*   Distance : " <<  distance[0] << G4endl;
      G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
             << G4endl;
#endif
          
   } else {  // p is on Hyperbolic surface.
   
      distance[0] = 0;
      xx.set(p.x(), p.y(), pz);
#ifdef __SOLIDDEBUG__
      G4cerr << "      ~~~ G4HyperbolicSurface::DistanceToSurface(p):~~~~~~" 
             << G4endl;
      G4cerr << "      //*   p is on of Hyperbolic surface."  << G4endl;
      G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
             << G4endl;
#endif
   }

   if (p.z() < 0) {
      G4ThreeVector tmpxx(xx.x(), xx.y(), -xx.z());
      xx = tmpxx;
   }

#ifdef __BOUNDARYCHECK__
   areacode[0] = GetAreaCode(xx, FALSE);
   if (!IsInside(areacode[0])) {
      // xx is out of boundary or corner.
      // return distance to boundary or corner.
      if (IsCorner(areacode[0])) {
         xx = GetCorner(areacode[0]);
         distance[0] = (xx - p).mag();
      } else {
         distance[0] = DistanceToBoundary(areacode[0], xx, p);
         G4cerr << "      G4HyperbolicSurface:DistanceToSurface(p) ~~~~~~~~~~~~~~" << G4endl;
         G4cerr << "         xx is out of boundary." << G4endl;
         G4cerr << "         areacode : " << std::hex << areacode[0] << std::dec << G4endl;
         G4cerr << "         xx : " << xx << G4endl;
         G4cerr << "         p : " << p << G4endl;
         G4cerr << "      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << G4endl;
         
      }
   }
#endif
       
   gxx[0] = ComputeGlobalPoint(xx);
   areacode[0]    = kInside;
   G4bool isvalid = TRUE;
   fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                             isvalid, 1, kDontValidate, &gp);
#ifdef __SOLIDDEBUG__
   G4cerr <<"      ~~~~~ G4HyperbolicSurface:DistanceToSurface(p):return"
          << G4endl;
   G4cerr <<"         NAME        : " << GetName() << G4endl;
   G4cerr <<"         xx, gxx     : " << xx << " " << gxx[0] << G4endl;
   G4cerr <<"         dist[0]     : " << distance[0] << G4endl;
   G4cerr <<"      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          << G4endl;
#endif

   return 1;
}

//=====================================================================
//* GetAreaCode -------------------------------------------------------
G4int G4HyperbolicSurface::GetAreaCode(const G4ThreeVector &xx, 
                                             G4bool         withTol)
{
   static const G4double ctol = 0.5 * kCarTolerance;
   G4int areacode = kInside;

   if ((fAxis[0] == kPhi && fAxis[1] == kZAxis))  {
      //G4int phiaxis = 0;
      G4int zaxis   = 1;
      
      if (withTol) {

         G4bool isoutside      = FALSE;
         G4int  phiareacode    = GetAreaCodeInPhi(xx);
         G4bool isoutsideinphi = IsOutside(phiareacode);

         // test boundary of phiaxis

         if ((phiareacode & kAxisMin) == kAxisMin) {

            areacode |= (kAxis0 & (kAxisPhi | kAxisMin)) | kBoundary;
            if (isoutsideinphi) isoutside = TRUE;

         } else if ((phiareacode & kAxisMax)  == kAxisMax) {

            areacode |= (kAxis0 & (kAxisPhi | kAxisMax)) | kBoundary;
            if (isoutsideinphi) isoutside = TRUE;

         }

         // test boundary of zaxis

         if (xx.z() < fAxisMin[zaxis] + ctol) {

            areacode |= (kAxis1 & (kAxisZ | kAxisMin));
            if   (areacode & kBoundary) areacode |= kCorner;  // xx is on the corner.
            else                        areacode |= kBoundary;

            if (xx.z() <= fAxisMin[zaxis] - ctol) isoutside = TRUE;

         } else if (xx.z() > fAxisMax[zaxis] - ctol) {

            areacode |= (kAxis1 & (kAxisZ | kAxisMax));
            if   (areacode & kBoundary) areacode |= kCorner;  // xx is on the corner.
            else                        areacode |= kBoundary;

            if (xx.z() >= fAxisMax[zaxis] + ctol) isoutside = TRUE;
         }

         // if isoutside = TRUE, clear kInside bit.
         // if not on boundary, add boundary information. 

         if (isoutside) {
            G4int tmpareacode = areacode & (~kInside);
            areacode = tmpareacode;
         } else if ((areacode & kBoundary) != kBoundary) {
            areacode |= (kAxis0 & kAxisPhi) | (kAxis1 & kAxisZ);
         }


#ifdef __SOLIDDEBUGAREACODE__
         G4cerr << "         === G4HyperbolicSurface::GetAreaCode ========="
                << G4endl;
         G4cerr << "         //#    NAME     : " << GetName() << G4endl;
         G4cerr << "         //#    xx       : " << xx << G4endl;
         G4cerr << "         //#    areacode : " << hex << areacode << G4endl;
         G4cerr << "         =============================================="
                << dec << G4endl;
#endif
         return areacode;
      
      } else {

         G4int phiareacode = GetAreaCodeInPhi(xx, FALSE);
         
         // test boundary of z-axis

         if (xx.z() < fAxisMin[zaxis]) {

            areacode |= (kAxis1 & (kAxisZ | kAxisMin)) | kBoundary;

         } else if (xx.z() > fAxisMax[zaxis]) {

            areacode |= (kAxis1 & (kAxisZ | kAxisMax)) | kBoundary;

         }

         // boundary of phi-axis

         if (phiareacode == kAxisMin) {

            areacode |= (kAxis0 & (kAxisPhi | kAxisMin));
            if   (areacode & kBoundary) areacode |= kCorner;  // xx is on the corner.
            else                        areacode |= kBoundary; 
             
         } else if (phiareacode == kAxisMax) {

            areacode |= (kAxis0 & (kAxisPhi | kAxisMax));
            if   (areacode & kBoundary) areacode |= kCorner;  // xx is on the corner.
            else                        areacode |= kBoundary; 
           
         }

         // if not on boundary, add boundary information. 

         if ((areacode & kBoundary) != kBoundary) {
            areacode |= (kAxis0 & kAxisPhi) | (kAxis1 & kAxisZ);
         }

         return areacode;

      }

   } else {
      G4cerr << "         G4HyperbolicSurface::GetAreaCode fAxis[0] = " 
             << fAxis[0] << " fAxis[1] = " << fAxis[1]
             << " is yet implemented. Write the code yourself." << G4endl;
      abort();
   }
}

//=====================================================================
//* GetAreaCodeInPhi --------------------------------------------------
G4int G4HyperbolicSurface::GetAreaCodeInPhi(const G4ThreeVector &xx,
                                                  G4bool withTol)
{
   
   G4ThreeVector lowerlimit; // lower phi-boundary limit at z = xx.z()
   G4ThreeVector upperlimit; // upper phi-boundary limit at z = xx.z()
   lowerlimit = GetBoundaryAtPZ(kAxis0 & kAxisMin, xx);
   upperlimit = GetBoundaryAtPZ(kAxis0 & kAxisMax, xx);

#ifdef __SOLIDDEBUGAREACODE__
   G4cerr << "         === G4HyperbolicSurface::GetAreaCodeInPhi ========="
          << G4endl;
   G4cerr << "         //#    NAME       : " << GetName() << G4endl;
   G4cerr << "         //#    xx         : " << xx << G4endl;
   G4cerr << "         //#    lowerlimit : " << lowerlimit << G4endl;
   G4cerr << "         //#    upperlimit : " << upperlimit << G4endl; 
   G4cerr << "         ==================================================="
          << G4endl;
 
#endif

   G4int  areacode  = kInside;
   G4bool isoutside = FALSE; 
   
   if (withTol) {
         
      if (AmIOnLeftSide(xx, lowerlimit) >= 0) {        // xx is on lowerlimit
         areacode |= (kAxisMin | kBoundary);
         if (AmIOnLeftSide(xx, lowerlimit) > 0) isoutside = TRUE; 

      } else if (AmIOnLeftSide(xx, upperlimit) <= 0) { // xx is on upperlimit
         areacode |= (kAxisMax | kBoundary);
         if (AmIOnLeftSide(xx, upperlimit) < 0) isoutside = TRUE; 
      }

      // if isoutside = TRUE, clear inside bit.

      if (isoutside) {
         G4int tmpareacode = areacode & (~kInside);
         areacode = tmpareacode;
      }


   } else {
   
      if (AmIOnLeftSide(xx, lowerlimit, FALSE) >= 0) {
         areacode |= (kAxisMin | kBoundary);
      } else if (AmIOnLeftSide(xx, upperlimit, FALSE) <= 0) {
         areacode |= (kAxisMax | kBoundary);
      }
   }

   return areacode;
   
}

//=====================================================================
//* SetCorners(solid) -------------------------------------------------
void G4HyperbolicSurface::SetCorners(G4TwistedTubs *solid)
{
   // Set Corner points in local coodinate.

   if (fAxis[0] == kPhi && fAxis[1] == kZAxis) {

      G4int i;
      G4double endPhi[2];
      G4double endRad[2];
      G4double endZ[2];
      G4double halfdphi = 0.5*(solid->GetDPhi());
      
      for (i=0; i<2; i++) { // i=0,1 : -ve z, +ve z
         endPhi[i] = solid->GetEndPhi(i);
         endZ[i]   = solid->GetEndZ(i);
         endRad[i] = (fHandedness == 1 ? solid->GetEndOuterRadius(i)
                                      : solid->GetEndInnerRadius(i));
      }

      G4int zmin = 0 ;  // at -ve z
      G4int zmax = 1 ;  // at +ve z

      G4double x, y, z;
      
      // corner of Axis0min and Axis1min
      x = endRad[zmin]*cos(endPhi[zmin] - halfdphi);
      y = endRad[zmin]*sin(endPhi[zmin] - halfdphi);
      z = endZ[zmin];
      SetCorner(kCorner0Min1Min, x, y, z);
      
      // corner of Axis0max and Axis1min
      x = endRad[zmin]*cos(endPhi[zmin] + halfdphi);
      y = endRad[zmin]*sin(endPhi[zmin] + halfdphi);
      z = endZ[zmin];
      SetCorner(kCorner0Max1Min, x, y, z);
      
      // corner of Axis0max and Axis1max
      x = endRad[zmax]*cos(endPhi[zmax] + halfdphi);
      y = endRad[zmax]*sin(endPhi[zmax] + halfdphi);
      z = endZ[zmax];
      SetCorner(kCorner0Max1Max, x, y, z);
      
      // corner of Axis0min and Axis1max
      x = endRad[zmax]*cos(endPhi[zmax] - halfdphi);
      y = endRad[zmax]*sin(endPhi[zmax] - halfdphi);
      z = endZ[zmax];
      SetCorner(kCorner0Min1Max, x, y, z);

   } else {
      G4cerr << "G4FlatSurface::SetCorners fAxis[0] = " << fAxis[0]
      << " fAxis[1] = " << fAxis[1]
      << " is yet implemented. Write the code yourself." << G4endl;
      abort();
   }
}

//=====================================================================
//* SetCorners() ------------------------------------------------------
void G4HyperbolicSurface::SetCorners()
{
   G4cerr << "G4FlatSurface::SetCorners"
   << " is yet implemented. Write the code yourself." << G4endl;
   abort();
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------
void G4HyperbolicSurface::SetBoundaries()
{
   // Set direction-unit vector of phi-boundary-lines in local coodinate.
   // kAxis0 must be kPhi.
   // This fanction set lower phi-boundary and upper phi-boundary.

   if (fAxis[0] == kPhi && fAxis[1] == kZAxis) {

      G4ThreeVector direction;
      // kAxis0 & kAxisMin
      direction = GetCorner(kCorner0Min1Max) - GetCorner(kCorner0Min1Min);
      direction = direction.unit();
      SetBoundary(kAxis0 & (kAxisPhi | kAxisMin), direction, 
                   GetCorner(kCorner0Min1Min), kAxisZ);

      // kAxis0 & kAxisMax
      direction = GetCorner(kCorner0Max1Max) - GetCorner(kCorner0Max1Min);
      direction = direction.unit();
      SetBoundary(kAxis0 & (kAxisPhi | kAxisMax), direction, 
                  GetCorner(kCorner0Max1Min), kAxisZ);

      // kAxis1 & kAxisMin
      direction = GetCorner(kCorner0Max1Min) - GetCorner(kCorner0Min1Min);
      direction = direction.unit();
      SetBoundary(kAxis1 & (kAxisZ | kAxisMin), direction, 
                   GetCorner(kCorner0Min1Min), kAxisPhi);

      // kAxis1 & kAxisMax
      direction = GetCorner(kCorner0Max1Max) - GetCorner(kCorner0Min1Max);
      direction = direction.unit();
      SetBoundary(kAxis1 & (kAxisZ | kAxisMax), direction, 
                  GetCorner(kCorner0Min1Max), kAxisPhi);
   } else {
      G4cerr << "G4HyperbolicSurface::SetBoundaries fAxis[0] = " << fAxis[0]
      << " fAxis[1] = " << fAxis[1]
      << " is yet implemented. Write the code yourself." << G4endl;
      abort();
   }
   
}





