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
// class G4ReplicaNavigation Implementation
//
// Author: P.Kent, 1996
//
// --------------------------------------------------------------------

#include "G4ReplicaNavigation.hh"

#include "G4AffineTransform.hh"
#include "G4SmartVoxelProxy.hh"
#include "G4SmartVoxelNode.hh"
#include "G4VSolid.hh"
#include "G4GeometryTolerance.hh"

namespace
{
  const G4ThreeVector VecCartAxes[3]=
  { G4ThreeVector(1.,0.,0.), G4ThreeVector(0.,1.,0.), G4ThreeVector(0.,0.,1.) };
  const G4ExitNormal::ESide SideCartAxesPlus[3]=
  { G4ExitNormal::kPX, G4ExitNormal::kPY, G4ExitNormal::kPZ };
  const G4ExitNormal::ESide SideCartAxesMinus[3]=
  { G4ExitNormal::kMX, G4ExitNormal::kMX, G4ExitNormal::kMX };
}

// ********************************************************************
// Constructor
// ********************************************************************
//
G4ReplicaNavigation::G4ReplicaNavigation()
{
  kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
  halfkCarTolerance = kCarTolerance*0.5;
  halfkRadTolerance = kRadTolerance*0.5;
  halfkAngTolerance = kAngTolerance*0.5;
  fMinStep = 0.05*kCarTolerance;
}

// ********************************************************************
// Destructor
// ********************************************************************
//
G4ReplicaNavigation::~G4ReplicaNavigation()
{
}

// ********************************************************************
// Inside
// ********************************************************************
//
EInside
G4ReplicaNavigation::Inside(const G4VPhysicalVolume* pVol,
                            const G4int replicaNo,
                            const G4ThreeVector& localPoint) const
{
  EInside in = kOutside;
  
  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;
  
  G4double coord, rad2, rmin, tolRMax2, rmax, tolRMin2;

  pVol->GetReplicationData(axis, nReplicas, width, offset, consuming);

  switch (axis)
  {
    case kXAxis:
    case kYAxis:
    case kZAxis:
      coord = std::fabs(localPoint(axis))-width*0.5;
      if ( coord<=-halfkCarTolerance )
      {
        in = kInside;
      }
      else if ( coord<=halfkCarTolerance )
      {
        in = kSurface;
      }
      break;
    case kPhi:
      if ( localPoint.y()||localPoint.x() )
      {
        coord = std::fabs(std::atan2(localPoint.y(),localPoint.x()))-width*0.5;
        if ( coord<=-halfkAngTolerance )
        {
          in = kInside;
        }
        else if ( coord<=halfkAngTolerance )
        {
          in = kSurface;
        }
      }
      else
      {
        in = kSurface;
      }
      break;
    case kRho:
      rad2 = localPoint.perp2();
      rmax = (replicaNo+1)*width+offset;
      tolRMax2  = rmax-halfkRadTolerance;
      tolRMax2 *= tolRMax2;
      if ( rad2>tolRMax2 )
      {
        tolRMax2 = rmax+halfkRadTolerance;
        tolRMax2 *= tolRMax2;
        if ( rad2<=tolRMax2 )
        {
          in = kSurface;
        }
      }
      else
      {
        // Known to be inside outer radius
        //
        if ( replicaNo||offset )
        {
          rmin = rmax-width;
          tolRMin2 = rmin-halfkRadTolerance;
          tolRMin2 *= tolRMin2;
          if ( rad2>tolRMin2 )
          {
            tolRMin2 = rmin+halfkRadTolerance;
            tolRMin2 *= tolRMin2;
            if ( rad2>=tolRMin2 )
            {
              in = kInside;
            }
            else
            {
              in = kSurface;
            }
          }
        }
        else
        {
          in = kInside;
        }
      }
      break;
    default:
      G4Exception("G4ReplicaNavigation::Inside()", "GeomNav0002",
                  FatalException, "Unknown axis!");
      break;
  }
  return in;
}

// ********************************************************************
// DistanceToOut
// ********************************************************************
//
G4double
G4ReplicaNavigation::DistanceToOut(const G4VPhysicalVolume* pVol,
                                   const G4int replicaNo,
                                   const G4ThreeVector& localPoint) const
{
  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width,offset;
  G4bool consuming;
  
  G4double safety = 0.;
  G4double safe1,safe2;
  G4double coord, rho, rmin, rmax;

  pVol->GetReplicationData(axis, nReplicas, width, offset, consuming);
  switch(axis)
  {
    case kXAxis:
    case kYAxis:
    case kZAxis:
       coord = localPoint(axis);
       safe1 = width*0.5-coord;
       safe2 = width*0.5+coord;
       safety = (safe1<=safe2) ? safe1 : safe2;
       break;
    case kPhi:
      if ( localPoint.y()<=0 )
      {
        safety = localPoint.x()*std::sin(width*0.5)
               + localPoint.y()*std::cos(width*0.5);
      }
      else
      {
        safety = localPoint.x()*std::sin(width*0.5)
               - localPoint.y()*std::cos(width*0.5);
      }
      break;
    case kRho:
      rho = localPoint.perp();
      rmax = width*(replicaNo+1)+offset;
      if ( replicaNo||offset )
      {
        rmin  = rmax-width;
        safe1 = rho-rmin;
        safe2 = rmax-rho;
        safety = (safe1<=safe2) ? safe1 : safe2;
      }
      else
      {
        safety = rmax-rho;
      }
      break;
    default:
     G4Exception("G4ReplicaNavigation::DistanceToOut()", "GeomNav0002",
                 FatalException, "Unknown axis!");
     break;
  }
  return (safety >= halfkCarTolerance) ? safety : 0;
}

// ********************************************************************
// DistanceToOut
// ********************************************************************
//
G4double
G4ReplicaNavigation::DistanceToOut(const G4VPhysicalVolume* pVol,
                                   const G4int replicaNo,
                                   const G4ThreeVector& localPoint,
                                   const G4ThreeVector& localDirection,
                                   G4ExitNormal& arExitNormal ) const
{
  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;

  G4double Dist=kInfinity;
  G4double coord, Comp, lindist;
  G4double signC = 0.0;
  G4ExitNormal candidateNormal; 
   
  pVol->GetReplicationData(axis, nReplicas, width, offset, consuming);
  switch(axis)
  {
    case kXAxis:
    case kYAxis:
    case kZAxis:
      coord = localPoint(axis);
      Comp = localDirection(axis);
      if ( Comp>0 )
      {
        lindist = width*0.5-coord;
        Dist = (lindist>0) ? lindist/Comp : 0;
        signC= 1.0;
      }
      else if ( Comp<0 )
      {
        lindist = width*0.5+coord;
        Dist = (lindist>0) ? -lindist/Comp : 0;
        signC= -1.0;
      }
      else
      {
        Dist = kInfinity;
      }
      // signC = sign<G4double>(Comp)
      candidateNormal.exitNormal = ( signC * VecCartAxes[axis]);
      candidateNormal.calculated = true;
      candidateNormal.validConvex = true;
      candidateNormal.exitSide =
        (Comp>0) ? SideCartAxesPlus[axis] : SideCartAxesMinus[axis];
      break;
    case kPhi:
      Dist = DistanceToOutPhi(localPoint,localDirection,width,candidateNormal);
        // candidateNormal set in call
      break;
    case kRho:
      Dist = DistanceToOutRad(localPoint,localDirection,width,offset,
                              replicaNo,candidateNormal);
        // candidateNormal set in call
      break;
    default:
     G4Exception("G4ReplicaNavigation::DistanceToOut()", "GeomNav0002",
                 FatalException, "Unknown axis!");
     break;
  }

  arExitNormal= candidateNormal; // .exitNormal;

  return Dist;
}

// ********************************************************************
// DistanceToOutPhi
// ********************************************************************
//
G4double
G4ReplicaNavigation::DistanceToOutPhi(const G4ThreeVector& localPoint,
                                      const G4ThreeVector& localDirection,
                                      const G4double width,
                                      G4ExitNormal& foundNormal ) const
{
  // Phi Intersection
  // NOTE: width<=pi by definition
  //
  G4double sinSPhi = -2.0, cosSPhi = -2.0;
  G4double pDistS, pDistE, compS, compE, Dist, dist2, yi;
  G4ExitNormal::ESide sidePhi = G4ExitNormal::kNull;
  G4ThreeVector  candidateNormal;

  if ( (localPoint.x()!=0.0) || (localPoint.y()!=0.0) )
  {
    sinSPhi = std::sin(-width*0.5);  // SIN of starting phi plane
    cosSPhi = std::cos(width*0.5);   // COS of starting phi plane

    // pDist -ve when inside
    //
    pDistS = localPoint.x()*sinSPhi-localPoint.y()*cosSPhi;
     // Start plane at phi= -S
    pDistE = localPoint.x()*sinSPhi+localPoint.y()*cosSPhi;
     // End   plane at phi= +S

    // Comp -ve when in direction of outwards normal
    //
    compS = -sinSPhi*localDirection.x()+cosSPhi*localDirection.y();
    compE = -sinSPhi*localDirection.x()-cosSPhi*localDirection.y();

    if ( (pDistS<=halfkCarTolerance)&&(pDistE<=halfkCarTolerance) )
    {
      // Inside both phi *full* planes
      //
      if ( compS<0 )
      {
        dist2 = pDistS/compS;
        yi = localPoint.y()+dist2*localDirection.y();

        // Check intersecting with correct half-plane (no -> no intersect)
        //
        if ( yi<=0 )
        {
          Dist = (pDistS<=-halfkCarTolerance) ? dist2 : 0;
          sidePhi= G4ExitNormal::kSPhi; // tbc
        }
        else
        {
          Dist = kInfinity;
        }
      }
      else
      {
        Dist = kInfinity;
      }
      if ( compE<0 )
      {
        dist2 = pDistE/compE;
        
        // Only check further if < starting phi intersection
        //
        if ( dist2<Dist )
        {
          yi = localPoint.y()+dist2*localDirection.y();

          // Check intersecting with correct half-plane
          //
          if ( yi>=0 )
          {
            // Leaving via ending phi
            //
            Dist = (pDistE<=-halfkCarTolerance) ? dist2 : 0;
            sidePhi = G4ExitNormal::kEPhi;
          }
        }
      }
    }
    else if ( (pDistS>halfkCarTolerance)&&(pDistE>halfkCarTolerance) )
    {
      // Outside both *full* phi planes
      // if towards both >=0 then once inside will remain inside
      //
      Dist = ((compS>=0)&&(compE>=0)) ? kInfinity : 0;
    }
    else if ( (pDistS>halfkCarTolerance)&&(pDistE<=halfkCarTolerance) )
    {
      // Outside full starting plane, inside full ending plane
      //
      if ( compE<0 )
      {      
        dist2 = pDistE/compE;
        yi = localPoint.y()+dist2*localDirection.y();

        // Check intersection in correct half-plane
        // (if not -> remain in extent)
        //
        Dist = (yi>0) ? dist2 : kInfinity;
        if( yi> 0 ) { sidePhi = G4ExitNormal::kEPhi; }
      }
      else  // Leaving immediately by starting phi
      {
        Dist = kInfinity;
      }
    }
    else
    {
      // Must be (pDistS<=halfkCarTolerance)&&(pDistE>halfkCarTolerance)
      // Inside full starting plane, outside full ending plane
      //
      if ( compE>=0 )
      {
        if ( compS<0 )
        {
          dist2 = pDistS/compS;
          yi = localPoint.y()+dist2*localDirection.y();

          // Check intersection in correct half-plane
          // (if not -> remain in extent)
          //
          Dist = (yi<0) ? dist2 : kInfinity;
          if(yi<0)  { sidePhi = G4ExitNormal::kSPhi; }
        }
        else
        {
          Dist = kInfinity;
        }
      }
      else
      {
        // Leaving immediately by ending phi
        //
        Dist = 0.;
        sidePhi = G4ExitNormal::kEPhi;
      }
    }
  }
  else
  {
    // On z axis + travel not || to z axis -> use direction vector
    //
    if( (std::fabs(localDirection.phi())<=width*0.5) )
    {
       Dist = kInfinity;
    }
    else
    {
       Dist = 0.;
       sidePhi = G4ExitNormal::kMY;
    }
  }

  if(sidePhi == G4ExitNormal::kSPhi )
  {
    candidateNormal = G4ThreeVector(sinSPhi,-cosSPhi,0.) ;
  }
  else if (sidePhi == G4ExitNormal::kEPhi)
  {
    candidateNormal = G4ThreeVector(sinSPhi,cosSPhi,0.) ;
  }
  else if (sidePhi == G4ExitNormal::kMY )
  {
    candidateNormal = G4ThreeVector(0., -1.0, 0.); // Split -S and +S 'phi'
  }
  foundNormal.calculated= (sidePhi != G4ExitNormal::kNull );
  foundNormal.exitNormal= candidateNormal;
   
  return Dist;
}

// ********************************************************************
// DistanceToOutRad
// ********************************************************************
//
G4double
G4ReplicaNavigation::DistanceToOutRad(const G4ThreeVector& localPoint,
                                      const G4ThreeVector& localDirection,
                                      const G4double width,
                                      const G4double offset,
                                      const G4int replicaNo,
                                      G4ExitNormal& foundNormal ) const
{
  G4double rmin, rmax, t1, t2, t3, deltaR;
  G4double b, c, d2, srd;
  G4ExitNormal::ESide  sideR= G4ExitNormal::kNull;

  //
  // Radial Intersections
  //
  
  // Find intersction with cylinders at rmax/rmin
  // Intersection point (xi,yi,zi) on line
  // x=localPoint.x+t*localDirection.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (localDirection.x^2+localDirection.y^2)t^2+
  //     2t(localPoint.x*localDirection.x+localPoint.y*localDirection.y)+
  //        localPoint.x^2+localPoint.y^2-R^2=0
  //
  //            t1                t2                    t3

  rmin = replicaNo*width+offset;
  rmax = (replicaNo+1)*width+offset;

  t1 = 1.0-localDirection.z()*localDirection.z();   // since v normalised
  t2 = localPoint.x()*localDirection.x()+localPoint.y()*localDirection.y();
  t3 = localPoint.x()*localPoint.x()+localPoint.y()*localPoint.y();
  
  if ( t1>0 )        // Check not parallel
  {
    // Calculate srd, r exit distance
    //
    if ( t2>=0 )
    {
      // Delta r not negative => leaving via rmax
      //
      deltaR = t3-rmax*rmax;
    
      // NOTE: Should use
      // rho-rmax<-halfkRadTolerance - [no sqrts for efficiency]
      //
      if ( deltaR<-halfkRadTolerance )
      {
        b  = t2/t1;
        c  = deltaR/t1;
        srd = -b+std::sqrt(b*b-c);
        sideR = G4ExitNormal::kRMax;
      }
      else
      {
        // On tolerant boundary & heading outwards (or locally
        // perpendicular to) outer radial surface -> leaving immediately
        //
        srd = 0;
        sideR = G4ExitNormal::kRMax;
      }
    }
    else
    {
      // Possible rmin intersection
      //
      if (rmin)
      {
        deltaR = t3-rmin*rmin;
        b  = t2/t1;
        c  = deltaR/t1;
        d2 = b*b-c;
        if ( d2>=0 )
        {
          // Leaving via rmin
          // NOTE: Should use
          // rho-rmin>halfkRadTolerance - [no sqrts for efficiency]
          //
          srd = (deltaR>halfkRadTolerance) ? -b-std::sqrt(d2) : 0.0;
          // Is the following more accurate ?
          // srd = (deltaR>halfkRadTolerance) ? c/( -b - std::sqrt(d2)) : 0.0;
          sideR = G4ExitNormal::kRMin;
        }
        else
        {
          // No rmin intersect -> must be rmax intersect
          //
          deltaR = t3-rmax*rmax;
          c  = deltaR/t1;
          d2 = b*b-c;
          srd = (d2 < 0.) ? 0.0 : -b+std::sqrt(d2);
          sideR = G4ExitNormal::kRMax;
        }
      }
      else
      {
        // No rmin intersect -> must be rmax intersect
        //
        deltaR = t3-rmax*rmax;
        b  = t2/t1;
        c  = deltaR/t1;
        d2 = b*b-c;
        srd = (d2 < 0.) ? 0.0 : -b+std::sqrt(d2);
        sideR= G4ExitNormal::kRMax;
      }
    }
  }
  else
  {
    srd =kInfinity;
    sideR = G4ExitNormal::kNull;
  }
   
  if( sideR != G4ExitNormal::kNull ) // if ((side == kRMax) || (side==kRMin))
  {
    // Note: returned vector not explicitly normalised
    // (divided by fRMax for unit vector)

    G4double xi, yi;
    xi = localPoint.x() + srd*localDirection.x();
    yi = localPoint.y() + srd*localDirection.y();
    G4ThreeVector normalR = G4ThreeVector(xi,yi,0.0);
     
    if( sideR == G4ExitNormal::kRMax )
    {
      normalR *= 1.0/rmax;
    }
    else
    {
      normalR *= (-1.0)/rmin;
    }
    foundNormal.exitNormal= normalR;
    foundNormal.calculated= true;
    foundNormal.validConvex = (sideR == G4ExitNormal::kRMax);
    foundNormal.exitSide = sideR;
  }
  else
  {
    foundNormal.calculated = false;
  }
   
  return srd;
}

// ********************************************************************
// ComputeTransformation
//
// Setup transformation and transform point into local system
// ********************************************************************
//
void
G4ReplicaNavigation::ComputeTransformation(const G4int replicaNo,
                                                 G4VPhysicalVolume* pVol,
                                                 G4ThreeVector& point) const
{
  G4double val,cosv,sinv,tmpx,tmpy;

  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width,offset;
  G4bool consuming;

  pVol->GetReplicationData(axis, nReplicas, width, offset, consuming);

  switch (axis)
  {
    case kXAxis:
      val = -width*0.5*(nReplicas-1)+width*replicaNo;
      pVol->SetTranslation(G4ThreeVector(val,0,0));
      point.setX(point.x()-val);
      break;
    case kYAxis:
      val = -width*0.5*(nReplicas-1)+width*replicaNo;
      pVol->SetTranslation(G4ThreeVector(0,val,0));
      point.setY(point.y()-val);
      break;
    case kZAxis:
      val = -width*0.5*(nReplicas-1)+width*replicaNo;
      pVol->SetTranslation(G4ThreeVector(0,0,val));
      point.setZ(point.z()-val);
      break;
    case kPhi:
      val = -(offset+width*(replicaNo+0.5));
      SetPhiTransformation(val,pVol);
      cosv = std::cos(val);
      sinv = std::sin(val);
      tmpx = point.x()*cosv-point.y()*sinv;
      tmpy = point.x()*sinv+point.y()*cosv;
      point.setY(tmpy);
      point.setX(tmpx);
      break;
    case kRho:
      // No setup required for radial case
    default:
      break;
  }
}

// ********************************************************************
// ComputeTransformation
//
// Setup transformation into local system
// ********************************************************************
//
void
G4ReplicaNavigation::ComputeTransformation(const G4int replicaNo,
                                                 G4VPhysicalVolume* pVol) const
{
  G4double val;

  // Replication data
  //
  EAxis axis;
  G4int nReplicas;
  G4double width, offset;
  G4bool consuming;

  pVol->GetReplicationData(axis, nReplicas, width, offset, consuming);

  switch (axis)
  {
    case kXAxis:
      val = -width*0.5*(nReplicas-1)+width*replicaNo;
      pVol->SetTranslation(G4ThreeVector(val,0,0));
      break;
    case kYAxis:
      val = -width*0.5*(nReplicas-1)+width*replicaNo;
      pVol->SetTranslation(G4ThreeVector(0,val,0));
      break;
    case kZAxis:
      val = -width*0.5*(nReplicas-1)+width*replicaNo;
      pVol->SetTranslation(G4ThreeVector(0,0,val));
      break;
    case kPhi:
      val = -(offset+width*(replicaNo+0.5));
      SetPhiTransformation(val,pVol);
      break;
    case kRho:
      // No setup required for radial case
    default:
      break;
  }
}

// ********************************************************************
// ComputeStep
// ********************************************************************
//
G4double
G4ReplicaNavigation::ComputeStep(const G4ThreeVector& globalPoint,
                                 const G4ThreeVector& globalDirection,
                                 const G4ThreeVector& localPoint,
                                 const G4ThreeVector& localDirection,
                                 const G4double currentProposedStepLength,
                                       G4double& newSafety,
                                       G4NavigationHistory &history,
                                 // std::pair<G4bool,G4bool>& validAndCalculated
                                       G4bool& validExitNormal,
                                       G4bool& calculatedExitNormal, 
                                       G4ThreeVector& exitNormalVector,
                                       G4bool& exiting,
                                       G4bool& entering,
                                       G4VPhysicalVolume* (*pBlockedPhysical),
                                       G4int& blockedReplicaNo )
{
  G4VPhysicalVolume *repPhysical, *motherPhysical;
  G4VPhysicalVolume *samplePhysical, *blockedExitedVol = nullptr;
  G4LogicalVolume *repLogical;
  G4VSolid *motherSolid;
  G4ThreeVector repPoint, repDirection, sampleDirection;
  G4double ourStep=currentProposedStepLength;
  G4double ourSafety=kInfinity;
  G4double sampleStep, sampleSafety, motherStep, motherSafety;
  G4long localNoDaughters, sampleNo;
  G4int depth;
  G4ExitNormal exitNormalStc;
  // G4int depthDeterminingStep= -1; // Useful only for debugging - for now

  calculatedExitNormal= false;
  
  // Exiting normal optimisation
  //
  if ( exiting&&validExitNormal )
  {
    if ( localDirection.dot(exitNormalVector)>=kMinExitingNormalCosine )
    {
      // Block exited daughter volume
      //
      blockedExitedVol = *pBlockedPhysical;
      ourSafety = 0;
    }
  }
  exiting  = false;
  entering = false;

  repPhysical = history.GetTopVolume();
  repLogical  = repPhysical->GetLogicalVolume();

  //
  // Compute intersection with replica boundaries & replica safety
  //

  sampleSafety = DistanceToOut(repPhysical,
                               history.GetTopReplicaNo(),
                               localPoint);
  G4ExitNormal normalOutStc;
  const G4int topDepth= (G4int)history.GetDepth();

  ourSafety = std::min( ourSafety, sampleSafety);

  if ( sampleSafety<ourStep )
  {

    sampleStep = DistanceToOut(repPhysical,
                               history.GetTopReplicaNo(),
                               localPoint,
                               localDirection,
                               normalOutStc);
    if ( sampleStep<ourStep )
    {
      ourStep = sampleStep;
      exiting = true;
      validExitNormal = normalOutStc.validConvex; // false; -> Old,Conservative

      exitNormalStc = normalOutStc;
      exitNormalStc.exitNormal =
        history.GetTopTransform().InverseTransformAxis(normalOutStc.exitNormal);
      calculatedExitNormal = true;
    }
  }
  const G4int secondDepth = topDepth;
  depth = secondDepth;  

  // Loop checking, 07.10.2016, JA -- Need to add: assert(depth>0)
  while ( history.GetVolumeType(depth)==kReplica )  
  {
    const G4AffineTransform& GlobalToLocal = history.GetTransform(depth);
    repPoint = GlobalToLocal.TransformPoint(globalPoint);
    // repPoint = history.GetTransform(depth).TransformPoint(globalPoint);
 
    sampleSafety = DistanceToOut(history.GetVolume(depth),
                                 history.GetReplicaNo(depth),
                                 repPoint);
    if ( sampleSafety < ourSafety )
    {
      ourSafety = sampleSafety;
    }
    if ( sampleSafety < ourStep )
    {
      G4ThreeVector newLocalDirection =
          GlobalToLocal.TransformAxis(globalDirection);
      sampleStep = DistanceToOut(history.GetVolume(depth),
                                 history.GetReplicaNo(depth),
                                 repPoint,
                                 newLocalDirection,
                                 normalOutStc);
      if ( sampleStep < ourStep )
      {
        ourStep = sampleStep;
        exiting = true;
       
        // As step is limited by this level, must set Exit Normal
        //
        G4ThreeVector localExitNorm = normalOutStc.exitNormal;
        G4ThreeVector globalExitNorm =
          GlobalToLocal.InverseTransformAxis(localExitNorm);

        exitNormalStc = normalOutStc; // Normal, convex, calculated, side
        exitNormalStc.exitNormal = globalExitNorm;
        calculatedExitNormal = true;
      }
    }
    depth--;
  }
 
  // Compute mother safety & intersection
  //
  G4ThreeVector exitVectorMother;
  G4bool exitConvex = false; // Value obtained in DistanceToOut(p,v) call
  G4ExitNormal motherNormalStc;

  repPoint = history.GetTransform(depth).TransformPoint(globalPoint);
  motherPhysical = history.GetVolume(depth);
  motherSolid = motherPhysical->GetLogicalVolume()->GetSolid();
  motherSafety = motherSolid->DistanceToOut(repPoint);
  repDirection = history.GetTransform(depth).TransformAxis(globalDirection);

  motherStep = motherSolid->DistanceToOut(repPoint,repDirection,true,
                                          &exitConvex,&exitVectorMother);
  if( exitConvex )
  {
     motherNormalStc = G4ExitNormal( exitVectorMother, true, false,
                                     G4ExitNormal::kMother);
     calculatedExitNormal = true;
  }
  const G4AffineTransform& globalToLocalTop = history.GetTopTransform();

  G4bool motherDeterminedStep = (motherStep<ourStep);

  if( (!exitConvex) && motherDeterminedStep )
  {
     exitVectorMother = motherSolid->SurfaceNormal( repPoint );
     motherNormalStc = G4ExitNormal( exitVectorMother, true, false,
                                     G4ExitNormal::kMother);
     // CalculatedExitNormal -> true;
     // Convex               -> false: do not know value
     // ExitSide             -> kMother (or kNull)
 
     calculatedExitNormal = true;
  }
  if( motherDeterminedStep )
  {
     G4ThreeVector globalExitNormalTop =
       globalToLocalTop.InverseTransformAxis(exitVectorMother);
     
     exitNormalStc = motherNormalStc;
     exitNormalStc.exitNormal = globalExitNormalTop;
  }

  // Push in principle no longer necessary. G4Navigator now takes care of ...
  // Removing this however may cause additional almost-zero steps and generate
  // warnings for pushed particles from G4Navigator, particularly for the case
  // of 3D replicas (Cartesian or combined Radial/Phi cases).
  // Requires further investigation and eventually reimplementation of
  // LevelLocate() to take into account point and direction ...
  //
  if( ourStep<fMinStep )
  {
    ourStep = 2*kCarTolerance;
  }

  if ( motherSafety<ourSafety )
  {
    ourSafety = motherSafety;
  }

#ifdef G4VERBOSE
  if ( fCheck )
  {
    if( motherSolid->Inside(localPoint)==kOutside )
    {
      std::ostringstream message;
      message << "Point outside volume !" << G4endl
              << "          Point " << localPoint
              << " is outside current volume " << motherPhysical->GetName()
              << G4endl;
      G4double estDistToSolid= motherSolid->DistanceToIn(localPoint); 
      message << "          Estimated isotropic distance to solid (distToIn)= " 
              << estDistToSolid << G4endl;
      if( estDistToSolid > 100.0 * kCarTolerance )
      {
        motherSolid->DumpInfo();
        G4Exception("G4ReplicaNavigation::ComputeStep()",
                    "GeomNav0003", FatalException, message,
                    "Point is far outside Current Volume !" ); 
      }
      else
        G4Exception("G4ReplicaNavigation::ComputeStep()",
                    "GeomNav1002", JustWarning, message,
                    "Point is a little outside Current Volume."); 
    }
  }
#endif

  // Comparison of steps may need precision protection
  //
#if 1
  if( motherDeterminedStep )
  {
    ourStep = motherStep;
    exiting = true;
  }

  // Transform it to the Grand-Mother Reference Frame (current convention)
  //
  if ( calculatedExitNormal )
  {
    if ( motherDeterminedStep )
    {
      exitNormalVector = motherNormalStc.exitNormal;
    }
    else
    {
      G4ThreeVector exitNormalGlobal = exitNormalStc.exitNormal;
      exitNormalVector = globalToLocalTop.TransformAxis(exitNormalGlobal);
      // exitNormalVector= globalToLocal2nd.TransformAxis(exitNormalGlobal);
      // Alt Make it in one go to Grand-Mother, avoiding transform below
    }
    // Transform to Grand-mother reference frame
    const G4RotationMatrix* rot = motherPhysical->GetRotation();
    if ( rot )
    {
      exitNormalVector *= rot->inverse();
    }

  }
  else
  {
    validExitNormal = false;
  }

#else
  if ( motherSafety<=ourStep )
  {
    if ( motherStep<=ourStep )
    {
      ourStep = motherStep;
      exiting = true;
      if ( validExitNormal )
      {
        const G4RotationMatrix* rot = motherPhysical->GetRotation();
        if ( rot )
        {
          exitNormal *= rot->inverse();
        }
      }
    }
    else
    {
      validExitNormal = false;
      // calculatedExitNormal= false;
    }
  }
#endif


  G4bool daughterDeterminedStep = false;
  G4ThreeVector daughtNormRepCrd;
     // Exit normal of daughter transformed to
     // the coordinate system of Replica (i.e. last depth)

  //
  // Compute daughter safeties & intersections
  //
  localNoDaughters = repLogical->GetNoDaughters();
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical = repLogical->GetDaughter((G4int)sampleNo);
    if ( samplePhysical!=blockedExitedVol )
    {
      G4ThreeVector localExitNorm;
      G4ThreeVector normReplicaCoord;

      G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                                 samplePhysical->GetTranslation());
      sampleTf.Invert();
      const G4ThreeVector samplePoint =
                        sampleTf.TransformPoint(localPoint);
      const G4VSolid* sampleSolid =
                        samplePhysical->GetLogicalVolume()->GetSolid();
      const G4double sampleSafetyDistance =
                        sampleSolid->DistanceToIn(samplePoint);
      if ( sampleSafetyDistance<ourSafety )
      {
        ourSafety = sampleSafetyDistance;
      }
      if ( sampleSafetyDistance<=ourStep )
      {
        sampleDirection = sampleTf.TransformAxis(localDirection);
        const G4double sampleStepDistance =
                        sampleSolid->DistanceToIn(samplePoint,sampleDirection);
        if ( sampleStepDistance<=ourStep )
        {
          daughterDeterminedStep = true;

          ourStep  = sampleStepDistance;
          entering = true;
          exiting  = false;
          *pBlockedPhysical = samplePhysical;
          blockedReplicaNo  = (G4int)sampleNo;

#ifdef DAUGHTER_NORMAL_ALSO
          // This norm can be calculated later, if needed daughter is available
          localExitNorm = sampleSolid->SurfaceNormal(samplePoint);
          daughtNormRepCrd = sampleTf.InverseTransformAxis(localExitNorm);
#endif
          
#ifdef G4VERBOSE
          // Check to see that the resulting point is indeed in/on volume.
          // This check could eventually be made only for successful candidate.

          if ( ( fCheck ) && ( sampleStepDistance < kInfinity ) )
          {
            G4ThreeVector intersectionPoint;
            intersectionPoint = samplePoint
                              + sampleStepDistance * sampleDirection;
            EInside insideIntPt = sampleSolid->Inside(intersectionPoint); 
            if ( insideIntPt != kSurface )
            {
              G4long oldcoutPrec = G4cout.precision(16); 
              std::ostringstream message;
              message << "Navigator gets conflicting response from Solid."
                      << G4endl
                      << "          Inaccurate DistanceToIn for solid "
                      << sampleSolid->GetName() << G4endl
                      << "          Solid gave DistanceToIn = "
                      << sampleStepDistance << " yet returns " ;
              if ( insideIntPt == kInside )
                message << "-kInside-"; 
              else if ( insideIntPt == kOutside )
                message << "-kOutside-";
              else
                message << "-kSurface-"; 
              message << " for this point !" << G4endl
                      << "          Point = " << intersectionPoint << G4endl;
              if ( insideIntPt != kInside )
                message << "        DistanceToIn(p) = " 
                       << sampleSolid->DistanceToIn(intersectionPoint)
                       << G4endl;
              if ( insideIntPt != kOutside ) 
                message << "        DistanceToOut(p) = " 
                       << sampleSolid->DistanceToOut(intersectionPoint);
              G4Exception("G4ReplicaNavigation::ComputeStep()", 
                          "GeomNav1002", JustWarning, message); 
              G4cout.precision(oldcoutPrec);
            }
          }
#endif
        }
      }
    }
  }

  calculatedExitNormal &= (!daughterDeterminedStep);

#ifdef DAUGHTER_NORMAL_ALSO
  if( daughterDeterminedStep )
  {
    // G4ThreeVector daughtNormGlobal =
    //   GlobalToLastDepth.Inverse().TransformAxis(daughtNormRepCrd);
    // ==> Can calculate it, but have no way to transmit it to caller (for now)

    exitNormalVector = globalToLocalTop.InverseTransformAxis(daughtNormGlobal);
    validExitNormal = false; // Entering daughter - never convex for parent

    calculatedExitNormal = true;
  }
  // calculatedExitNormal= true;  // Force it to true -- dubious
#endif

  newSafety = ourSafety;
  return ourStep;
}

// ********************************************************************
// ComputeSafety
//
// Compute the isotropic distance to current volume's boundaries
// and to daughter volumes.
// ********************************************************************
//
G4double
G4ReplicaNavigation::ComputeSafety(const G4ThreeVector& globalPoint,
                                   const G4ThreeVector& localPoint,
                                         G4NavigationHistory& history,
                                   const G4double )
{
  G4VPhysicalVolume *repPhysical, *motherPhysical;
  G4VPhysicalVolume *samplePhysical, *blockedExitedVol = nullptr;
  G4LogicalVolume *repLogical;
  G4VSolid *motherSolid;
  G4ThreeVector repPoint;
  G4double ourSafety = kInfinity;
  G4double sampleSafety;
  G4long localNoDaughters, sampleNo;
  G4int depth;

  repPhysical = history.GetTopVolume();
  repLogical  = repPhysical->GetLogicalVolume();

  //
  // Compute intersection with replica boundaries & replica safety
  //

  sampleSafety = DistanceToOut(history.GetTopVolume(),
                               history.GetTopReplicaNo(),
                               localPoint);
  if ( sampleSafety<ourSafety )
  {
    ourSafety = sampleSafety;
  }

  depth = (G4int)history.GetDepth()-1;

  // Loop checking, 07.10.2016, JA -- need to add: assert(depth>0)
  while ( history.GetVolumeType(depth)==kReplica )
  {      
    repPoint = history.GetTransform(depth).TransformPoint(globalPoint);
    sampleSafety = DistanceToOut(history.GetVolume(depth),
                                 history.GetReplicaNo(depth),
                                 repPoint);
    if ( sampleSafety<ourSafety )
    {
      ourSafety = sampleSafety;
    }
    depth--;
  }

  // Compute mother safety & intersection
  //
  repPoint = history.GetTransform(depth).TransformPoint(globalPoint);
  motherPhysical = history.GetVolume(depth);
  motherSolid = motherPhysical->GetLogicalVolume()->GetSolid();
  sampleSafety = motherSolid->DistanceToOut(repPoint);

  if ( sampleSafety<ourSafety )
  {
    ourSafety = sampleSafety;
  }

  // Compute daughter safeties & intersections
  //
  localNoDaughters = repLogical->GetNoDaughters();
  for ( sampleNo=localNoDaughters-1; sampleNo>=0; sampleNo-- )
  {
    samplePhysical = repLogical->GetDaughter((G4int)sampleNo);
    if ( samplePhysical!=blockedExitedVol )
    {
      G4AffineTransform sampleTf(samplePhysical->GetRotation(),
                                 samplePhysical->GetTranslation());
      sampleTf.Invert();
      const G4ThreeVector samplePoint =
                            sampleTf.TransformPoint(localPoint);
      const G4VSolid *sampleSolid =
                            samplePhysical->GetLogicalVolume()->GetSolid();
      const G4double sampleSafetyDistance =
                            sampleSolid->DistanceToIn(samplePoint);
      if ( sampleSafetyDistance<ourSafety )
      {
        ourSafety = sampleSafetyDistance;
      }
    }
  }
  return ourSafety;
}

// ********************************************************************
// BackLocate
// ********************************************************************
//
EInside
G4ReplicaNavigation::BackLocate(G4NavigationHistory& history,
                          const G4ThreeVector& globalPoint,
                                G4ThreeVector& localPoint,
                          const G4bool& exiting,
                                G4bool& notKnownInside ) const
{
  G4VPhysicalVolume *pNRMother = nullptr;
  G4VSolid *motherSolid;
  G4ThreeVector repPoint, goodPoint;
  G4int mdepth, depth, cdepth;
  EInside insideCode;

  cdepth = (G4int)history.GetDepth();
  
  // Find non replicated mother
  //
  for ( mdepth=cdepth-1; mdepth>=0; mdepth-- )
  {
    if ( history.GetVolumeType(mdepth)!=kReplica )
    {
      pNRMother = history.GetVolume(mdepth);
      break;
    }
  }

  if( pNRMother == nullptr ) 
  {
    // All the tree of mother volumes were Replicas. 
    // This is an error, as the World volume must be a Placement
    //
    G4Exception("G4ReplicaNavigation::BackLocate()", "GeomNav0002",
                FatalException, "The World volume must be a Placement!");
    return kInside;
  }

  motherSolid = pNRMother->GetLogicalVolume()->GetSolid();
  goodPoint = history.GetTransform(mdepth).TransformPoint(globalPoint);
  insideCode = motherSolid->Inside(goodPoint);
  if ( (insideCode==kOutside)||((insideCode==kSurface)&&exiting) )
  {
    // Outside mother -> back up to mother level
    // Locate.. in Navigator will back up one more level
    // localPoint not required
    //
    history.BackLevel(cdepth-mdepth);
    //      localPoint = goodPoint;
  }
  else
  {
    notKnownInside = false;

    // Still within replications
    // Check down: if on outside stop at this level
    //
    for ( depth=mdepth+1; depth<cdepth; ++depth)
    {
      repPoint = history.GetTransform(depth).TransformPoint(globalPoint);
      insideCode = Inside(history.GetVolume(depth),
                          history.GetReplicaNo(depth),
                          repPoint);
      if ( (insideCode==kOutside)||((insideCode==kSurface)&&exiting) )
      {
        localPoint = goodPoint;
        history.BackLevel(cdepth-depth);
        return insideCode;
      }
      else
      {
        goodPoint = repPoint;
      }
    }
    localPoint = history.GetTransform(depth).TransformPoint(globalPoint);
    insideCode = Inside(history.GetVolume(depth),
                        history.GetReplicaNo(depth),
                        localPoint);
    // If outside level, set localPoint = coordinates in reference system
    // of *previous* level - location code in navigator will back up one
    // level [And also manage blocking]
    //
    if ( (insideCode==kOutside)||((insideCode==kSurface)&&exiting) )
    {
      localPoint = goodPoint;
    }
  }
  return insideCode;
}
