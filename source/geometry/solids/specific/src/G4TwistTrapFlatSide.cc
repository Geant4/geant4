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
// $Id: G4TwistTrapFlatSide.cc 66356 2012-12-18 09:02:32Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4TwistTrapFlatSide.cc
//
// Author: 
//   30-Aug-2002 - O.Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include "G4TwistTrapFlatSide.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4TwistTrapFlatSide::G4TwistTrapFlatSide( const G4String        &name,
                              G4double      PhiTwist,
                              G4double      pDx1,
                              G4double      pDx2,
                              G4double      pDy,
                              G4double      pDz,
                              G4double      pAlpha,
                              G4double      pPhi,
                              G4double      pTheta,
                              G4int         handedness) 

  : G4VTwistSurface(name)
{
   fHandedness = handedness;   // +z = +ve, -z = -ve

   fDx1 = pDx1 ;
   fDx2 = pDx2 ;
   fDy = pDy ;
   fDz = pDz ;
   fAlpha = pAlpha ;
   fTAlph = std::tan(fAlpha) ;
   fPhi  = pPhi ;
   fTheta = pTheta ;

   fdeltaX = 2 * fDz * std::tan(fTheta) * std::cos(fPhi)  ;
     // dx in surface equation
   fdeltaY = 2 * fDz * std::tan(fTheta) * std::sin(fPhi)  ;
     // dy in surface equation

   fPhiTwist = PhiTwist ;

   fCurrentNormal.normal.set( 0, 0, (fHandedness < 0 ? -1 : 1)); 
         // Unit vector, in local coordinate system
   fRot.rotateZ( fHandedness > 0 
                 ? 0.5 * fPhiTwist
                 : -0.5 * fPhiTwist );

   fTrans.set(
              fHandedness > 0 ? 0.5*fdeltaX : -0.5*fdeltaX , 
              fHandedness > 0 ? 0.5*fdeltaY : -0.5*fdeltaY ,
              fHandedness > 0 ? fDz : -fDz ) ;

   fIsValidNorm = true;


   fAxis[0] = kXAxis ;
   fAxis[1] = kYAxis ;
   fAxisMin[0] = kInfinity ;  // x-Axis cannot be fixed, because it 
   fAxisMax[0] =  kInfinity ; // depends on y
   fAxisMin[1] = -fDy ;  // y - axis
   fAxisMax[1] =  fDy ;

   SetCorners();
   SetBoundaries();
}


//=====================================================================
//* Fake default constructor ------------------------------------------

G4TwistTrapFlatSide::G4TwistTrapFlatSide( __void__& a )
  : G4VTwistSurface(a), fDx1(0.), fDx2(0.), fDy(0.), fDz(0.), fPhiTwist(0.), 
    fAlpha(0.), fTAlph(0.), fPhi(0.), fTheta(0.), fdeltaX(0.), fdeltaY(0.)
{
}


//=====================================================================
//* destructor --------------------------------------------------------

G4TwistTrapFlatSide::~G4TwistTrapFlatSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4TwistTrapFlatSide::GetNormal(const G4ThreeVector & /* xx */ , 
                                             G4bool isGlobal)
{
   if (isGlobal) {
      return ComputeGlobalDirection(fCurrentNormal.normal);
   } else {
      return fCurrentNormal.normal;
   }
}

//=====================================================================
//* DistanceToSurface(p, v) -------------------------------------------

G4int G4TwistTrapFlatSide::DistanceToSurface(const G4ThreeVector &gp,
                                       const G4ThreeVector &gv,
                                             G4ThreeVector  gxx[],
                                             G4double       distance[],
                                             G4int          areacode[],
                                             G4bool         isvalid[],
                                             EValidate      validate) 
{
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

   //
   // special case!
   // if p is on surface, distance = 0. 
   //

   if (std::fabs(p.z()) == 0.) {   // if p is on the plane
      distance[0] = 0;
      G4ThreeVector xx = p;
      gxx[0] = ComputeGlobalPoint(xx);
      
      if (validate == kValidateWithTol) {
         areacode[0] = GetAreaCode(xx);
         if (!IsOutside(areacode[0])) {
            isvalid[0] = true;
         }
      } else if (validate == kValidateWithoutTol) {
         areacode[0] = GetAreaCode(xx, false);
         if (IsInside(areacode[0])) {
            isvalid[0] = true;
         }
      } else { // kDontValidate
         areacode[0] = sInside;
         isvalid[0] = true;
      }

      return 1;
   }
   //
   // special case end
   //
   
   if (v.z() == 0) { 

      fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0], 
                                     isvalid[0], 0, validate, &gp, &gv);
      return 0;
   }
   
   distance[0] = - (p.z() / v.z());
      
   G4ThreeVector xx = p + distance[0]*v;
   gxx[0] = ComputeGlobalPoint(xx);

   if (validate == kValidateWithTol) {
      areacode[0] = GetAreaCode(xx);
      if (!IsOutside(areacode[0])) {
         if (distance[0] >= 0) isvalid[0] = true;
      }
   } else if (validate == kValidateWithoutTol) {
      areacode[0] = GetAreaCode(xx, false);
      if (IsInside(areacode[0])) {
         if (distance[0] >= 0) isvalid[0] = true;
      }
   } else { // kDontValidate
      areacode[0] = sInside;
         if (distance[0] >= 0) isvalid[0] = true;
   }


   fCurStatWithV.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                                  isvalid[0], 1, validate, &gp, &gv);

#ifdef G4TWISTDEBUG
   G4cerr << "ERROR - G4TwistTrapFlatSide::DistanceToSurface(p,v)" << G4endl;
   G4cerr << "        Name        : " << GetName() << G4endl;
   G4cerr << "        xx          : " << xx << G4endl;
   G4cerr << "        gxx[0]      : " << gxx[0] << G4endl;
   G4cerr << "        dist[0]     : " << distance[0] << G4endl;
   G4cerr << "        areacode[0] : " << areacode[0] << G4endl;
   G4cerr << "        isvalid[0]  : " << isvalid[0]  << G4endl;
#endif
   return 1;
}

//=====================================================================
//* DistanceToSurface(p) ----------------------------------------------

G4int G4TwistTrapFlatSide::DistanceToSurface(const G4ThreeVector &gp,
                                             G4ThreeVector  gxx[],
                                             G4double       distance[],
                                             G4int          areacode[])
{
   // Calculate distance to plane in local coordinate,
   // then return distance and global intersection points.
   //  

   fCurStat.ResetfDone(kDontValidate, &gp);

   if (fCurStat.IsDone()) {
      G4int i;
      for (i=0; i<fCurStat.GetNXX(); i++) {
         gxx[i] = fCurStat.GetXX(i);
         distance[i] = fCurStat.GetDistance(i);
         areacode[i] = fCurStat.GetAreacode(i);
      }
      return fCurStat.GetNXX();
   } else {
      // initialize
      G4int i;
      for (i=0; i<2; i++) {
         distance[i] = kInfinity;
         areacode[i] = sOutside;
         gxx[i].set(kInfinity, kInfinity, kInfinity);
      }
   }
   
   G4ThreeVector p = ComputeLocalPoint(gp);
   G4ThreeVector xx;

   // The plane is placed on origin with making its normal 
   // parallel to z-axis. 
   if (std::fabs(p.z()) <= 0.5 * kCarTolerance)
   {   // if p is on the plane, return 1
      distance[0] = 0;
      xx = p;
   } else {
      distance[0] = std::fabs(p.z());
      xx.set(p.x(), p.y(), 0);  
   }

   gxx[0] = ComputeGlobalPoint(xx);
   areacode[0] = sInside;
   G4bool isvalid = true;
   fCurStat.SetCurrentStatus(0, gxx[0], distance[0], areacode[0],
                             isvalid, 1, kDontValidate, &gp);
   return 1;

}

G4int G4TwistTrapFlatSide::GetAreaCode(const G4ThreeVector &xx, 
                                       G4bool withTol)
{

  static const G4double ctol = 0.5 * kCarTolerance;
  G4int areacode = sInside;
  
  if (fAxis[0] == kXAxis && fAxis[1] == kYAxis) {

    G4int yaxis = 1;
    
    G4double wmax = xAxisMax(xx.y(), fTAlph ) ;
    G4double wmin = -xAxisMax(xx.y(), -fTAlph ) ;

    if (withTol) {
      
      G4bool isoutside   = false;
      
      // test boundary of x-axis
      
      if (xx.x() < wmin + ctol) {
        areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary; 
        if (xx.x() <= wmin - ctol) isoutside = true;
        
      } else if (xx.x() > wmax - ctol) {
        areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
        if (xx.x() >= wmax + ctol)  isoutside = true;
      }
      
      // test boundary of y-axis
      
      if (xx.y() < fAxisMin[yaxis] + ctol) {
        areacode |= (sAxis1 & (sAxisY | sAxisMin)); 
        
        if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
        else                        areacode |= sBoundary;
        if (xx.y() <= fAxisMin[yaxis] - ctol) isoutside = true;
        
      } else if (xx.y() > fAxisMax[yaxis] - ctol) {
        areacode |= (sAxis1 & (sAxisY | sAxisMax));
        
        if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
        else                        areacode |= sBoundary; 
        if (xx.y() >= fAxisMax[yaxis] + ctol) isoutside = true;
      }
      
      // if isoutside = true, clear inside bit.             
      // if not on boundary, add axis information.             
      
      if (isoutside) {
        G4int tmpareacode = areacode & (~sInside);
        areacode = tmpareacode;
      } else if ((areacode & sBoundary) != sBoundary) {
        areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisY);
      }           
      
    } else {
      
      // boundary of x-axis
      
      if (xx.x() < wmin ) {
        areacode |= (sAxis0 & (sAxisX | sAxisMin)) | sBoundary;
      } else if (xx.x() > wmax) {
        areacode |= (sAxis0 & (sAxisX | sAxisMax)) | sBoundary;
      }
      
      // boundary of y-axis
      
      if (xx.y() < fAxisMin[yaxis]) {
        areacode |= (sAxis1 & (sAxisY | sAxisMin));
        if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
        else                        areacode |= sBoundary; 
        
      } else if (xx.y() > fAxisMax[yaxis]) {
        areacode |= (sAxis1 & (sAxisY | sAxisMax)) ;
        if   (areacode & sBoundary) areacode |= sCorner;  // xx is on the corner.
        else                        areacode |= sBoundary; 
      }
      
      if ((areacode & sBoundary) != sBoundary) {
        areacode |= (sAxis0 & sAxisX) | (sAxis1 & sAxisY);
      }           
    }
    return areacode;
  } else {
    G4Exception("G4TwistTrapFlatSide::GetAreaCode()",
                "GeomSolids0001", FatalException,
                "Feature NOT implemented !");
  }
  
  return areacode;
}


//=====================================================================
//* SetCorners --------------------------------------------------------

void G4TwistTrapFlatSide::SetCorners()
{
   // Set Corner points in local coodinate.

   if (fAxis[0] == kXAxis && fAxis[1] == kYAxis) {
   
     G4double x, y, z;
      
     // corner of Axis0min and Axis1min
     x = -fDx1 + fDy * fTAlph ;
     y = -fDy ;
     z = 0 ;
     SetCorner(sC0Min1Min, x, y, z);
      
     // corner of Axis0max and Axis1min
     x = fDx1 + fDy * fTAlph ;
     y = -fDy ;
     z = 0 ;
     SetCorner(sC0Max1Min, x, y, z);
     
     // corner of Axis0max and Axis1max
     x = fDx2 - fDy * fTAlph ;
     y = fDy ;
     z = 0 ;
     SetCorner(sC0Max1Max, x, y, z);
     
     // corner of Axis0min and Axis1max
     x = -fDx2 - fDy * fTAlph ;
     y = fDy ;
     z = 0 ;
     SetCorner(sC0Min1Max, x, y, z);
     
   } else {
     std::ostringstream message;
     message << "Feature NOT implemented !" << G4endl
             << "        fAxis[0] = " << fAxis[0] << G4endl
             << "        fAxis[1] = " << fAxis[1];
     G4Exception("G4TwistTrapFlatSide::SetCorners()",
                 "GeomSolids0001", FatalException, message);
   }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4TwistTrapFlatSide::SetBoundaries()
{
   // Set direction-unit vector of phi-boundary-lines in local coodinate.
   // Don't call the function twice.
   
  G4ThreeVector direction ;

  if (fAxis[0] == kXAxis && fAxis[1] == kYAxis) {
   
    // sAxis0 & sAxisMin
    direction = - ( GetCorner(sC0Min1Max) - GetCorner(sC0Min1Min) ) ;
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisX | sAxisMin), direction, 
                GetCorner(sC0Min1Max), sAxisY) ;
    
    // sAxis0 & sAxisMax
    direction = GetCorner(sC0Max1Max) - GetCorner(sC0Max1Min)  ; // inverse
    direction = direction.unit();
    SetBoundary(sAxis0 & (sAxisX | sAxisMax), direction, 
                GetCorner(sC0Max1Min), sAxisY);
    
    // sAxis1 & sAxisMin
    direction = GetCorner(sC0Max1Min) - GetCorner(sC0Min1Min);
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisY | sAxisMin), direction, 
                GetCorner(sC0Min1Min), sAxisX);
    
    // sAxis1 & sAxisMax
    direction = - ( GetCorner(sC0Max1Max) - GetCorner(sC0Min1Max) ) ;
    direction = direction.unit();
    SetBoundary(sAxis1 & (sAxisY | sAxisMax), direction, 
                GetCorner(sC0Max1Max), sAxisX);
    
  } else {
    std::ostringstream message;
    message << "Feature NOT implemented !" << G4endl
            << "        fAxis[0] = " << fAxis[0] << G4endl
            << "        fAxis[1] = " << fAxis[1];
    G4Exception("G4TwistTrapFlatSide::SetCorners()",
                "GeomSolids0001", FatalException, message);
  }
}

//=====================================================================
//* GetFacets() -------------------------------------------------------

void G4TwistTrapFlatSide::GetFacets( G4int k, G4int n, G4double xyz[][3],
                                     G4int faces[][4], G4int iside ) 
{

  G4double x,y    ;     // the two parameters for the surface equation
  G4ThreeVector p ;  // a point on the surface, given by (z,u)

  G4int nnode ;
  G4int nface ;

  G4double xmin,xmax ;

  // calculate the (n-1)*(k-1) vertices

  G4int i,j ;

  for ( i = 0 ; i<n ; i++ ) {

    y = -fDy + i*(2*fDy)/(n-1) ;

    for ( j = 0 ; j<k ; j++ ) {

      xmin = GetBoundaryMin(y) ;
      xmax = GetBoundaryMax(y) ;
      x = xmin + j*(xmax-xmin)/(k-1) ;

      nnode = GetNode(i,j,k,n,iside) ;
      p = SurfacePoint(x,y,true) ;  // surface point in global coordinate system

      xyz[nnode][0] = p.x() ;
      xyz[nnode][1] = p.y() ;
      xyz[nnode][2] = p.z() ;

      if ( i<n-1 && j<k-1 ) {   

        nface = GetFace(i,j,k,n,iside) ;

        if (fHandedness < 0) {  // lower side 
          faces[nface][0] = GetEdgeVisibility(i,j,k,n,0,1) * ( GetNode(i  ,j  ,k,n,iside)+1) ;  
          faces[nface][1] = GetEdgeVisibility(i,j,k,n,1,1) * ( GetNode(i+1,j  ,k,n,iside)+1) ;
          faces[nface][2] = GetEdgeVisibility(i,j,k,n,2,1) * ( GetNode(i+1,j+1,k,n,iside)+1) ;
          faces[nface][3] = GetEdgeVisibility(i,j,k,n,3,1) * ( GetNode(i  ,j+1,k,n,iside)+1) ;
        } else {                // upper side
          faces[nface][0] = GetEdgeVisibility(i,j,k,n,0,-1) * ( GetNode(i  ,j  ,k,n,iside)+1) ;  
          faces[nface][1] = GetEdgeVisibility(i,j,k,n,1,-1) * ( GetNode(i  ,j+1,k,n,iside)+1) ;
          faces[nface][2] = GetEdgeVisibility(i,j,k,n,2,-1) * ( GetNode(i+1,j+1,k,n,iside)+1) ;
          faces[nface][3] = GetEdgeVisibility(i,j,k,n,3,-1) * ( GetNode(i+1,j  ,k,n,iside)+1) ;
        }

      }
    }
  }
}
