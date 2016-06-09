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
// $Id: G4FlatTrapSide.cc,v 1.5 2004/12/08 10:20:37 link Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4FlatTrapSide.cc
//
// Author: 
//   30-Aug-2002 - Oliver Link (Oliver.Link@cern.ch)
//
// --------------------------------------------------------------------

#include "G4FlatTrapSide.hh"

//=====================================================================
//* constructors ------------------------------------------------------

G4FlatTrapSide::G4FlatTrapSide( const G4String        &name,
                              G4double      PhiTwist,
                              G4double      pDx1,
                              G4double      pDx2,
                              G4double      pDy,
                              G4double      pDz,
                              G4int         handedness) 

  : G4VSurface(name)
{
   fHandedness = handedness;   // +z = +ve, -z = -ve

   fDx1 = pDx1 ;
   fDx2 = pDx2 ;
   fDy = pDy ;
   fDz = pDz ;

   fPhiTwist = PhiTwist ;

   fCurrentNormal.normal.set(0, 0, (fHandedness < 0 ? -1 : 1)); 
         // Unit vector, in local coordinate system
   fRot.rotateZ( fHandedness > 0 
                 ? 0.5 * fPhiTwist
                 : -0.5 * fPhiTwist );

   fTrans.set(0, 0, fHandedness > 0 ? fDz : -fDz ) ;

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
//* destructor --------------------------------------------------------

G4FlatTrapSide::~G4FlatTrapSide()
{
}

//=====================================================================
//* GetNormal ---------------------------------------------------------

G4ThreeVector G4FlatTrapSide::GetNormal(const G4ThreeVector & /* xx */ , 
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

G4int G4FlatTrapSide::DistanceToSurface(const G4ThreeVector &gp,
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

#ifdef G4SPECSDEBUG
   G4cerr << "ERROR - G4FlatTrapSide::DistanceToSurface(p,v)" << G4endl;
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

G4int G4FlatTrapSide::DistanceToSurface(const G4ThreeVector &gp,
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
   if (std::fabs(p.z()) <= 0.5 * kCarTolerance) {   // if p is on the plane, return 1
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

G4int G4FlatTrapSide::GetAreaCode(const G4ThreeVector &xx, 
                                       G4bool withTol)
{

  static const G4double ctol = 0.5 * kCarTolerance;
  G4int areacode = sInside;
  
  if (fAxis[0] == kXAxis && fAxis[1] == kYAxis) {

    G4int yaxis = 1;
    
   G4double wmax = fDx2 + ( fDx1-fDx2)/2. - xx.y() * (fDx1-fDx2)/(2*fDy) ;
   G4double wmin = -wmax ;


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
    G4Exception("G4FlatTrapSide::GetAreaCode()",
                "NotImplemented", FatalException,
                "Feature NOT implemented !");
  }
  
  return areacode;
}


//=====================================================================
//* SetCorners --------------------------------------------------------

void G4FlatTrapSide::SetCorners()
{
   // Set Corner points in local coodinate.

   if (fAxis[0] == kXAxis && fAxis[1] == kYAxis) {
   
     G4double x, y, z;
      
     // corner of Axis0min and Axis1min
     x = -fDx1 ;
     y = -fDy ;
     z = 0 ;
     SetCorner(sC0Min1Min, x, y, z);
      
     // corner of Axis0max and Axis1min
     x = fDx1 ;
     y = -fDy ;
     z = 0 ;
     SetCorner(sC0Max1Min, x, y, z);
     
     // corner of Axis0max and Axis1max
     x = fDx2 ;
     y = fDy ;
     z = 0 ;
     SetCorner(sC0Max1Max, x, y, z);
     
     // corner of Axis0min and Axis1max
     x = -fDx2 ;
     y = fDy ;
     z = 0 ;
     SetCorner(sC0Min1Max, x, y, z);
     
   } else {
     G4cerr << "ERROR - G4FlatTrapSide::SetCorners()" << G4endl
            << "        fAxis[0] = " << fAxis[0] << G4endl
            << "        fAxis[1] = " << fAxis[1] << G4endl;
     G4Exception("G4FlatTrapSide::SetCorners()",
                 "NotImplemented", FatalException,
                 "Feature NOT implemented !");
   }
}

//=====================================================================
//* SetBoundaries() ---------------------------------------------------

void G4FlatTrapSide::SetBoundaries()
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
    G4cerr << "ERROR - G4FlatTrapSide::SetBoundaries()" << G4endl
           << "        fAxis[0] = " << fAxis[0] << G4endl
           << "        fAxis[1] = " << fAxis[1] << G4endl;
    G4Exception("G4FlatTrapSide::SetCorners()",
                "NotImplemented", FatalException,
                "Feature NOT implemented !");
   }
}
