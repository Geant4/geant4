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

// $Id: testG4Sphere2.cc,v 1.1 2003-11-13 19:33:31 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4Sphere Test File
//
// o Basic asserts on each function +
//   awkward cases for tracking / geom algorithms
//
// o Add tests on dicovering bugs in G4Sphere.cc...
//
// History:
// 28.03.95 P.Kent Initial version
// 20.10.96 V.Grichine Final modifications to commit

#include "G4ios.hh"
#include <assert.h>
#include <math.h>
#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4Sphere.hh"

//const G4double kApproxEqualTolerance = kCarTolerance;

// Return true if the double check is approximately equal to target
//
// Process:
//
// Return true is difference < kApproxEqualTolerance

//G4bool ApproxEqual(const G4double check,const G4double target)
//{
//    return (fabs(check-target)<kApproxEqualTolerance) ? true : false ;
//}

// Return true if the 3vector check is approximately equal to target
//G4bool ApproxEqual(const G4ThreeVector& check, const G4ThreeVector& target)
//{
//    return (ApproxEqual(check.x(),target.x())&&
//	   ApproxEqual(check.y(),target.y())&&
//	    ApproxEqual(check.z(),target.z()))? true : false;
//}

///////////////////////////////////////////////////////////////////
//
// Dave's auxiliary function

const G4String OutputInside(const EInside a)
{
	switch(a) 
        {
		case kInside:  return "Inside"; 
		case kOutside: return "Outside";
		case kSurface: return "Surface";
	}
	return "????";
}



int main(void)
{
    G4double Dist;
    G4ThreeVector pzero(0,0,0),px(30,0,0),py(0,30,0),pz(0,0,30);
    G4ThreeVector Pmx(-30,0,0),pmy(0,-30,0),pmz(0,0,-30);
    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector ponrmin1(45,0,0),ponrmax1(50,0,0),ponzmax(0,0,50),
	    ponrmin2(45/sqrt(2.),45/sqrt(2.),0),
            ponrmin3(0,0,-45),ponrminJ(0,0,-300),ponrmaxJ(0,0,-500),
	    ponrmax2(50/sqrt(2.),50/sqrt(2.),0);
    G4ThreeVector ponphi1(48/sqrt(2.),-48/sqrt(2.),0),
	          ponphi2(48/sqrt(2.),48/sqrt(2.),0),
	          pInPhi(48*0.866,-24,0),
	          pOverPhi(-48/sqrt(2.),48/sqrt(2.),0);
    G4ThreeVector pontheta1(0,48*sin(pi/4),48*cos(pi/4)),
	    pontheta2(0,48*sin(pi/4),-48*cos(pi/4));

    G4ThreeVector ptestphi1(-100,-45/sqrt(2.),0),
	    ptestphi2(-100,45/sqrt(2.),0);

    G4ThreeVector ptesttheta1(0,48/sqrt(2.),100),
	    ptesttheta2(0,48/sqrt(2.),-100);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/sqrt(2.),1/sqrt(2.),0),vmxmy(-1/sqrt(2.),-1/sqrt(2.),0);
    G4ThreeVector vxmy(1/sqrt(2.),-1/sqrt(2.),0),vmxy(-1/sqrt(2.),1/sqrt(2.),0);
    G4ThreeVector v345exit1(-0.8,0.6,0),v345exit2(0.8,0.6,0),
	          v345exit3(0.6,0.8,0);
    G4ThreeVector norm,*pNorm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;


        G4bool checkPoint( const G4Sphere& pSph, G4ThreeVector origin,
                           G4double  d,    G4ThreeVector dir,    EInside  exp); 

        G4Sphere SpAroundX("SpAroundX",  10.*mm, 1000.*mm, -1.0*degree, 2.0*degree, 0.*degree, 180.0*degree );

	G4double  sinOneDeg = sin( 1.0 * degree );
	G4double  radOne = 100.0 * mm;

	G4ThreeVector  ptPhiSurfExct= G4ThreeVector( radOne * cos( -1.0 * degree ) , 
			             - radOne *  sinOneDeg, 
				      0.0 );
        G4cout << " Starting from point " << ptPhiSurfExct << G4endl;
        G4cout << "   using direction   " << vy << G4endl; 
        checkPoint( SpAroundX, ptPhiSurfExct, -radOne * kAngTolerance * 1.5, 
                    vy,   kOutside); 
        checkPoint( SpAroundX, ptPhiSurfExct, -radOne * kAngTolerance * 0.49, 
                    vy,   kSurface); 
        checkPoint( SpAroundX, ptPhiSurfExct, -radOne * kAngTolerance * 0.25, 
                    vy,   kSurface); 
        checkPoint( SpAroundX, ptPhiSurfExct,  0.0,
                    vy,   kSurface); 
        checkPoint( SpAroundX, ptPhiSurfExct,  radOne * kAngTolerance * 0.25, 
                    vy,   kSurface); 
        checkPoint( SpAroundX, ptPhiSurfExct,  radOne * kAngTolerance * 0.49, 
                    vy,   kSurface); 
        checkPoint( SpAroundX, ptPhiSurfExct,  radOne * kAngTolerance * 1.5, 
                    vy,   kInside); 

        // Try one that has a 'deep' negative phi section
	//   --> Vlad. Grichine test case, 30 Oct 2003
        // 
        G4cout << G4endl << G4endl << "" << G4endl;
	G4cout << "========================================================= " << G4endl; 

	G4Sphere SphDeepNeg("DeepNegPhiSphere",  10.*mm, 1000.*mm, 
			   -270.0*degree, 280.0*degree,          //  start Phi,   delta Phi
			    0.*degree, 180.0*degree );        //  start Theta, delta Theta
        G4double phiPoint = 160.0 * degree; 
        G4ThreeVector  StartPt( radOne * cos(phiPoint), radOne * sin(phiPoint), 0.0); 
        G4cout << "For sphere " << SphDeepNeg.GetName() << G4endl;
        G4cout << " Starting from point " << ptPhiSurfExct << G4endl;

        checkPoint( SphDeepNeg, StartPt,  0.0,  vy,   kInside); 

        // Try the edges  
        G4ThreeVector  NegEdgePt( radOne * cos(-270.0*degree), radOne * sin(-270.0*degree), 0.0); 
        G4ThreeVector  PosEdgePt( radOne * cos(10.0*degree), radOne * sin(10.0*degree), 0.0); 

        G4cout << "--------------------------------------------------------" << G4endl; 
	G4cout << " New point " << NegEdgePt << " should be at Neg edge of -270.0 degrees " <<  G4endl;
	checkPoint( SphDeepNeg, NegEdgePt,  0.0,  -vx,   kSurface); 
	checkPoint( SphDeepNeg, NegEdgePt,  radOne*kAngTolerance * 0.25,  -vx,   kSurface); 
	checkPoint( SphDeepNeg, NegEdgePt, -radOne*kAngTolerance * 0.25,  -vx,   kSurface); 
	checkPoint( SphDeepNeg, NegEdgePt,  radOne*kAngTolerance * 1.25,  -vx,   kInside); 
	checkPoint( SphDeepNeg, NegEdgePt, -radOne*kAngTolerance * 1.25,  -vx,   kOutside); 

	G4cout << "--------------------------------------------------------" << G4endl; 
	G4cout << " New point " << PosEdgePt << " should be at Pos edge of +10.0 degrees " <<  G4endl;
        checkPoint( SphDeepNeg, PosEdgePt,  0.0,  -vy,   kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt,  radOne*kAngTolerance * 0.25,  -vy,   kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt, -radOne*kAngTolerance * 0.25,  -vy,   kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt, -radOne*kAngTolerance * 1.25,  -vy,   kOutside); 
	checkPoint( SphDeepNeg, PosEdgePt,  radOne*kAngTolerance * 1.25,  -vy,   kInside); 

        G4double radMax= 1000.0 * mm; 
	NegEdgePt = G4ThreeVector( radMax * cos(-270.0*degree), radMax * sin(-270.0*degree), 0.0); 
        G4cout << "--------------------------------------------------------" << G4endl; 
	G4cout << " New point " << NegEdgePt << " should be at RadMax / Neg edge of -270.0 degrees " <<  G4endl;
	checkPoint( SphDeepNeg, NegEdgePt,  0.0,  -vx,   kSurface); 
	checkPoint( SphDeepNeg, NegEdgePt,  radMax*kAngTolerance * 0.25,  -vx,   kSurface); 
	checkPoint( SphDeepNeg, NegEdgePt, -radMax*kAngTolerance * 0.25,  -vx,   kSurface); 
	checkPoint( SphDeepNeg, NegEdgePt,  radMax*kAngTolerance * 1.25,  -vx,   kSurface); 
	checkPoint( SphDeepNeg, NegEdgePt, -radMax*kAngTolerance * 1.25,  -vx,   kOutside); 

        PosEdgePt= G4ThreeVector( radMax * cos(10.0*degree), radMax * sin(10.0*degree), 0.0); 
	G4cout << "--------------------------------------------------------" << G4endl; 
	G4cout << " New point " << PosEdgePt << " should be at RadMax Pos edge of +10.0 degrees " <<  G4endl;
        checkPoint( SphDeepNeg, PosEdgePt,  0.0,  -vy,   kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt,  radMax*kAngTolerance * 0.25,  -vy,   kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt, -radMax*kAngTolerance * 0.25,  -vy,   kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt, -radMax*kAngTolerance * 1.25,  -vy,   kOutside); 
	checkPoint( SphDeepNeg, PosEdgePt,  radMax*kAngTolerance * 1.25,  -vy,   kSurface); 

	G4ThreeVector NormInDir = - cos(10.0*degree) * vy + sin(10.0*degree) * vx; 
        checkPoint( SphDeepNeg, PosEdgePt,  0.0,  -vy,   kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt,  radMax*kAngTolerance * 0.25, NormInDir, kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt, -radMax*kAngTolerance * 0.25, NormInDir, kSurface); 
	checkPoint( SphDeepNeg, PosEdgePt, -radMax*kAngTolerance * 1.25, NormInDir, kOutside); 
	checkPoint( SphDeepNeg, PosEdgePt,  radMax*kAngTolerance * 1.25, NormInDir, kSurface); 

	return 0;
}

// Given the sphere 'rSphere', a point 'origin'
//                                               -->     -->              --> 
//   the function below checks that the point    pnew=( origin + dist * direction )
//     - the point (p+dist*dir) is located in the part of the solid given by 'expectedInResult'
//     - and that from there, the DistanceToIn along 'dir' is not negative or Infinite
//
//  Use cases expected:
//   - 'origin' is on/near a surface and 'direction' is pointing towards the inside of the solid

G4bool
checkPoint( const G4Sphere &rSphere, 
            G4ThreeVector  origin, 
            G4double       dist, 
            G4ThreeVector  direction, 
            EInside        expectedInResult)
{
    G4int verbose = 0, verboseErr= 2; 
    G4bool  testAll = false ;   //  if false, do not call DistToIn if Outside etc.

    G4ThreeVector newPoint; 
    G4double  distIn=-1.0, distOut=-1.0; 

    newPoint = origin + dist * direction; 

    G4int oldPrecision= G4cout.precision(10); 
    // G4cout << " --- Sphere " << rSphere.GetName() << "" << G4endl;
    if( verbose > 0 ) {
      G4cout << G4endl;
      if (verbose > 2 ) G4cout << " Sphere " << rSphere.GetName();
      G4cout.precision(10); 
      if (verbose > 1 ) G4cout << " dir= " << direction;
      G4cout << " dist= " << dist;
    } 

    EInside  inSphere=  rSphere.Inside( newPoint ) ; 
                              /*======*/
    G4cout.precision(15); 
    // G4cout << " NewPoint  " << newPoint << " is " 
    G4bool goodIn= (inSphere == expectedInResult) ; 
    if ( !goodIn ) {
      G4cout << " ************ Unexpected Result for Inside *************** " << G4endl;
    } 
    if ( verbose || !goodIn ) {
       G4cout << " New point " 
	      << " is "  <<  OutputInside( inSphere ) 
	      <<  " vs " <<  OutputInside( expectedInResult ) 
	      <<  " expected." <<  G4endl ;
    }

    G4bool goodDistIn = true; 

    distIn = rSphere.DistanceToIn( newPoint, direction ); 
                    /*===========*/
    if ( verbose )  G4cout << " DistToIn (p, dir) = " << distIn << G4endl;
    if( (inSphere == kOutside) 
        && (distIn < 0.0 ) // Cannot use 0.5*kCarTolerance for Angular tolerance!! 
      ){
       G4cout << " ********** Unexpected Result for DistanceToIn from outside ********* " << G4endl;
       // G4cout << " It should be " << G4endl;
       goodDistIn = false;
    }
    if( (inSphere == kSurface ) 
        && ( (distIn < 0.0) || (distIn >= kInfinity )) 
      ){
       G4cout << " ********** Unexpected Result for DistanceToIn on surface   ********* " << G4endl;
       // if ( (distIn != 0.0) ) 
       //  -  Can check that the return value must be 0.0
       //     But in general case the direction can be away from the solid, 
       //        and then a finite or kInfinity answer is correct
       //        --> must check the direction against the normal
       //            in order to perform this check in general case.

       goodDistIn = false;
    }
    if ( verbose || !goodDistIn ) {
      G4cout << " DistToIn (p, dir) = " << distIn << G4endl;
    }

    G4bool good= (goodIn && goodDistIn);  
    if ( !good ){ 
       // There was an error -- document the use case!
       G4cout << " --- Sphere " << rSphere.GetName() << "" << G4endl;
       if ( verboseErr > 1 ) { 
	 G4cout << "  dist= " << dist ; // << G4endl;
	 G4cout << "  Direction= " << direction  ; // << G4endl; 
	 G4cout << "  dist/kAngTolerance = " << dist / kAngTolerance << G4endl;
  	 G4cout << "  Original pt=  " << origin << G4endl; 
     }
       G4cout <<   "  Actual-point= " << newPoint << G4endl;
       G4cout <<   "  Rho= " << G4ThreeVector(newPoint.x(), newPoint.y(), 0.).mag() << G4endl;     
    } 

    if( testAll || (inSphere!=kOutside) ) { 
      
       distOut = rSphere.DistanceToOut( newPoint, direction ); 
                       /*=============*/
       if ( verbose ) G4cout << " DistToOut (p, dir) = " << distOut << G4endl;
    }

    G4cout.precision(oldPrecision); 

    return good;
}
