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

// $Id$
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
#undef NDEBUG
#include <assert.h>
#include <cmath>
#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "ApproxEqual.hh"
#include "G4GeometryTolerance.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4Sphere.hh"
#include "Randomize.hh"

//const G4double kApproxEqualTolerance = kCarTolerance;

// Return true if the double check is approximately equal to target
//
// Process:
//
// Return true is difference < kApproxEqualTolerance

//G4bool ApproxEqual(const G4double check,const G4double target)
//{
//    return (std::fabs(check-target)<kApproxEqualTolerance) ? true : false ;
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
    G4double kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();
    EInside inside;
    G4int i, iMax;
    G4double Dist, vol, volCheck;
    G4double phi, cosTheta, sinTheta, rMax, rRand, zP;

    G4ThreeVector pzero(0,0,0),px(30,0,0),py(0,30,0),pz(0,0,30);
    G4ThreeVector pmx(-30,0,0),pmy(0,-30,0),pmz(0,0,-30);
    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector ponrmin1(45,0,0),ponrmax1(50,0,0),ponzmax(0,0,50),
	    ponrmin2(45/std::sqrt(2.),45/std::sqrt(2.),0),
            ponrmin3(0,0,-45),ponrminJ(0,0,-300),ponrmaxJ(0,0,-500),
	    ponrmax2(50/std::sqrt(2.),50/std::sqrt(2.),0);
    G4ThreeVector ponphi1(48/std::sqrt(2.),-48/std::sqrt(2.),0),
	          ponphi2(48/std::sqrt(2.),48/std::sqrt(2.),0),
	          pInPhi(48*0.866,-24,0),
	          pOverPhi(-48/std::sqrt(2.),48/std::sqrt(2.),0);
    G4ThreeVector pontheta1(0,48*std::sin(pi/4),48*std::cos(pi/4)),
	    pontheta2(0,48*std::sin(pi/4),-48*std::cos(pi/4));

    G4ThreeVector ptestphi1(-100,-45/std::sqrt(2.),0),
	    ptestphi2(-100,45/std::sqrt(2.),0);

    G4ThreeVector ptesttheta1(0,48/std::sqrt(2.),100),
	    ptesttheta2(0,48/std::sqrt(2.),-100);

    // Directions

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/std::sqrt(2.),1/std::sqrt(2.),0),vmxmy(-1/std::sqrt(2.),-1/std::sqrt(2.),0);
    G4ThreeVector vxmy(1/std::sqrt(2.),-1/std::sqrt(2.),0),vmxy(-1/std::sqrt(2.),1/std::sqrt(2.),0);
    G4ThreeVector vxmz(1/std::sqrt(2.),0.,-1/std::sqrt(2.)),vymz(0.,1/std::sqrt(2.),-1/std::sqrt(2.));
    G4ThreeVector vmxmz(-1/std::sqrt(2.),0.,-1/std::sqrt(2.)),vmxz(-1/std::sqrt(2.),0.,1/std::sqrt(2.));

    G4ThreeVector v345exit1(-0.8,0.6,0),v345exit2(0.8,0.6,0),
	          v345exit3(0.6,0.8,0);





    G4ThreeVector pRand, vRand, norm, *pNorm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4Sphere s1("Solid G4Sphere",0,50,0,twopi,0,pi);
    G4Sphere sn1("sn1",0,50,halfpi,3.*halfpi,0,pi);

    // Theta sections

    G4Sphere sn11("sn11",0,50,0,twopi,0.,halfpi);
    G4Sphere sn12("sn12",0,50,0,twopi,0.,0.25*pi);
    G4Sphere sn13("sn13",0,50,0,twopi,0.75*pi,0.25*pi);
    G4Sphere sn14("sn14",0,50,0,twopi,0.25*pi,0.75*pi);
    G4Sphere sn15("sn15",0,50,0,twopi,89.*deg,91.*deg);
    G4Sphere sn16("sn16",0,50,0,twopi,0.*deg,89.*deg);
    G4Sphere sn17("sn17",0,50,0,twopi,91.*deg,89.*deg);

    G4Sphere s2("Spherical Shell",45,50,0,twopi,0,pi);
    G4Sphere sn2("sn2",45,50,halfpi,halfpi,0,pi);
    G4Sphere sn22("sn22",0,50,halfpi,halfpi,0,pi);



    G4Sphere s3("Band (theta segment)",45,50,0,twopi,pi/4,halfpi);
    G4Sphere s32("Band (theta segment2)",45,50,0,twopi,0,pi/4);
    G4Sphere s33("Band (theta segment1)",45,50,0,twopi,pi*3/4,pi/4);
    G4Sphere s34("Band (theta segment)",4,50,0,twopi,pi/4,halfpi);
    G4Sphere s4("Band (phi segment)",45,50,-pi/4,halfpi,0,twopi);
    //    G4cout<<"s4.fSPhi = "<<s4.GetSPhi()<<G4endl;
    G4Sphere s41("Band (phi segment)",5,50,-pi,3.*pi/2.,0,twopi);
    G4Sphere s42("Band (phi segment)",5,50,-pi/2,3.*pi/2.,0,twopi);
    G4Sphere s5("Patch (phi/theta seg)",45,50,-pi/4,halfpi,pi/4,halfpi);

    G4Sphere s6("John example",300,500,0,5.76,0,pi) ; 
    G4Sphere s7("sphere7",1400.,1550.,0.022321428571428572,0.014642857142857141,
	                  1.5631177553663251,0.014642857142857141    );
    G4Sphere s8("sphere",278.746*mm, 280.0*mm, 0.0*degree, 360.0*degree,
		                               0.0*degree, 90.0*degree);

    G4Sphere b216("b216", 1400.0, 1550.0, 
                  0.022321428571428572, 
		  0.014642857142857141,
                  1.578117755366325,
                  0.014642857142867141);



    G4Sphere s9("s9",0*mm,410*mm,0*degree,360*degree,90*degree,90*degree);

    G4Sphere b402("b402", 475*mm, 480*mm, 
               0*degree,360*degree,17.8*degree,144.4*degree);

    G4Sphere b1046("b1046", 4750*km, 4800*km, 
               0*degree,360*degree,0*degree,90*degree);

    G4Sphere testTheta("opticalEscape",0.,2.995,0*degree,360*degree,0*degree,120*degree);
   
G4ThreeVector p402(471.7356120367253*mm, 51.95081450791341*mm, 5.938043020529463*mm);
G4ThreeVector v402(-0.519985502840818, 0.2521089719986221, 0.8161226274728446);



G4ThreeVector p216(1549.9518578505142,1.2195415370970153,-12.155289555510985), 
              v216(-0.61254821852534425,-0.51164551429243466,-0.60249775741147549);


G4ThreeVector s9p(384.8213314370455*mm,
	       134.264386151667*mm,
	       -44.56026800002064*mm);

G4ThreeVector s9v(-0.6542770611918751,
		   -0.0695116921641141,
		   -0.7530535517814154);

  G4Sphere s10("s10",0*mm,0.018*mm,0*degree,360*degree,0*degree,180*degree);

  G4ThreeVector s10p(0.01160957408065766*mm,
                     0.01308205826682229*mm,0.004293345210644617*mm);

  G4Sphere s11("s11",5000000.*mm,
                    3700000000.*mm,
                   0*degree,360*degree,0*degree,180*degree);

 // G4ThreeVector ps11(-1184000000.*mm,-3477212676.843337059*mm,-444000000.*mm);

 // G4ThreeVector ps11( -3339032195.112830162*mm, -1480000000*mm, -592000000*mm );

  G4ThreeVector ps11( -3072559844.81995153427124*mm, -1924000000*mm, -740000000*mm );

  G4Sphere sAlex("sAlex",500.*mm,
                    501.*mm,
                   0*degree,360*degree,0*degree,180*degree);

  G4ThreeVector psAlex(-360.4617031263808*mm,
                     -158.1198807105035*mm,308.326878333183*mm);

  G4ThreeVector vsAlex(-0.7360912456240805,-0.4955800202572754,0.4610532741813497 );


  G4Sphere sLHCB("sLHCB",8600*mm, 8606*mm, 
    -1.699135525184141*degree,
    3.398271050368263*degree,88.52855940538514*degree,2.942881189229715*degree );

  G4ThreeVector pLHCB(8600.242072535835,-255.1193517702246,-69.0010277128286);



  G4Sphere b658("b658", 209.6*mm, 211.2658*mm, 
                           0.0*degree, 360*degree, 0.0*degree, 90*degree); 


  G4ThreeVector p658(-35.69953348982516*mm, 198.3183279249958, 56.30959457033987 ); 
  G4ThreeVector v658(-.2346058124516908,-0.9450502890785083,0.2276841318065671); 


  G4Sphere spAroundX("SpAroundX",  10.*mm, 1000.*mm, -1.0*degree, 
                                                     2.0*degree, 
                                        0.*degree, 180.0*degree );

  G4double  radOne = 100.0*mm;
  G4double  angle = -1.0*degree - 0.25*kAngTolerance;
  G4ThreeVector  ptPhiMinus= G4ThreeVector( radOne*std::cos(angle) ,
                                           radOne*std::sin(angle),
                                           0.0 );

  // spheres for theta cone intersections

  G4Sphere s13("Band (theta segment)",5,50,0,twopi,pi/6.,halfpi);
  G4Sphere s14("Band (theta segment)",5,50,0,twopi,pi/3.,halfpi);


  // b. 830

  G4double mainInnerRadius =  21.45 * cm ;
  G4double mainOuterRadius = 85.0 * cm ;
  G4double minTheta = 18.0 * degree ;
 

  
 
   G4Sphere sb830( "mainSp",
                    mainInnerRadius,
                    mainOuterRadius,
                    0.0, M_PI*2,
                    minTheta,
                    M_PI - 2*minTheta);

   
G4ThreeVector pb830(81.61117212,-27.77179755,196.4143423);
   G4ThreeVector vb830(0.1644697995,0.18507236,0.9688642354);


   G4Sphere      s18_03_10("s18_03_10", 47*mm, 48*mm, 0*deg, 200*deg, 80*deg, 100*deg);
   G4ThreeVector p18_03_10(43.37539710867407*mm, 16.12036885157033*mm, -8.224548565698871*mm);
   G4ThreeVector v18_03_10(-0.4161958548132512, 0.6603942936714858, -0.6250283092167956);

   inside = s18_03_10.Inside(p18_03_10) ;
   G4cout<<"s18_03_10.Inside(p18_03_10 ... = "<<OutputInside(inside)<<G4endl ;



#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif

   G4cout.precision(20);

    //////////////// Check name /////////////////////////

    assert(s1.GetName()=="Solid G4Sphere");
    
    // check cubic volume

    vol = s1.GetCubicVolume();
    volCheck = 4*pi*125000/3; 

    assert( ApproxEqual(vol,volCheck));

    assert(ApproxEqual(s1.GetCubicVolume(),4*pi*125000/3));


    // Some user application cases

    Dist = s4.DistanceToOut(ponrmin3,vmz,calcNorm,pgoodNorm,pNorm) ;
    // G4cout<<"s4.DistanceToOut(ponrmin3,vmz,calcNorm,pgoodNorm,pNorm) = "<<Dist
    //    <<G4endl ;
    Dist = s6.DistanceToOut(ponrminJ,vmz,calcNorm,pgoodNorm,pNorm) ;
    // G4cout<<"s6.DistanceToOut(ponrminJ,vmz,calcNorm,pgoodNorm,pNorm) = "<<Dist
    //    <<G4endl ;
    Dist = s6.DistanceToOut(ponrmaxJ,vmz,calcNorm,pgoodNorm,pNorm) ;
    // G4cout<<"s6.DistanceToOut(ponrmaxJ,vmz,calcNorm,pgoodNorm,pNorm) = "<<Dist
    //  <<G4endl ;

    Dist = s7.DistanceToOut(G4ThreeVector(1399.984667238032,
                                             5.9396696802500299,
                                            -2.7661927818688308),
                            G4ThreeVector(0.50965504781062942,
                                         -0.80958849145715217,
                                          0.29123565499656401),
                                        calcNorm,pgoodNorm,pNorm) ;
    // G4cout<<"s7.DistanceToOut(shereP,sphereV,calcNorm,pgoodNorm,pNorm) = "<<Dist
    //   <<G4endl ;

    Dist = s7.DistanceToIn(G4ThreeVector(1399.984667238032,
                                             5.9396696802500299,
                                            -2.7661927818688308),
                            G4ThreeVector(0.50965504781062942,
                                         -0.80958849145715217,
                                          0.29123565499656401)) ;
    // G4cout<<"s7.DistanceToIn(shereP,sphereV,calcNorm,pgoodNorm,pNorm) = "<<Dist
    //    <<G4endl ;


    
// Check G4Sphere::Inside

    inside = s7.Inside(G4ThreeVector(1399.984667238032,
                                             5.9396696802500299,
                                            -2.7661927818688308)  ) ;
    // G4cout<<"s7.Inside(G4ThreeVector(1399.98466 ... = "
    //       <<OutputInside(inside)<<G4endl ;

    inside = s8.Inside(G4ThreeVector(-249.5020724528353*mm,
					       26.81253142743162*mm,
                                              114.8988524453591*mm  )  ) ;
    // G4cout<<"s8.Inside(G4ThreeVector(-249.5020 ... = "<<OutputInside(inside)<<G4endl ;
    inside = b216.Inside(p216);
    // G4cout<<"b216.Inside(p216) = "<<OutputInside(inside)<<G4endl ;
    inside = s1.Inside(pz);
    // G4cout<<"s1.Inside(pz) = "<<OutputInside(inside)<<G4endl ;

    inside = s9.Inside(s9p);
    // G4cout<<"s9.Inside(s9p) = "<<OutputInside(inside)<<G4endl ;

    inside = b402.Inside(p402);
    // G4cout<<"p402.Inside(p402) = "<<OutputInside(inside)<<G4endl ;

    inside = s10.Inside(s10p);
    // G4cout<<"s10.Inside(s10p) = "<<OutputInside(inside)<<G4endl ;
    // G4cout<<"p radius = "<<s10p.mag()<<G4endl ;
    
    inside = s11.Inside(ps11);
    // G4cout<<"s11.Inside(ps11) = "<<OutputInside(inside)<<G4endl ;
    // G4cout<<"ps11.mag() = "<<ps11.mag()<<G4endl ;
    
    inside = sLHCB.Inside(pLHCB);
    // G4cout<<"sLHCB.Inside(pLHCB) = "<<OutputInside(inside)<<G4endl ;
    // G4cout<<"pLHCB.mag() = "<<pLHCB.mag()<<G4endl ;
    
    inside = spAroundX.Inside(ptPhiMinus);
    // G4cout<<"spAroundX.Inside(ptPhiMinus) = "<<OutputInside(inside)<<G4endl ;
     inside = b658.Inside(p658);
    G4cout<<"b658.Inside(p658) = "<<OutputInside(inside)<<G4endl ;

    assert(s1.Inside(pzero)==kInside);
    // assert(s1.Inside(pz)==kInside);
    assert(s2.Inside(pzero)==kOutside);
    assert(s2.Inside(ponrmin2)==kSurface);
    assert(s2.Inside(ponrmax2)==kSurface);
    assert(s3.Inside(pontheta1)==kSurface);
    assert(s3.Inside(pontheta2)==kSurface);
    assert(s4.Inside(ponphi1)==kSurface);
    assert(s4.Inside(ponphi1)==kSurface);
    assert(s4.Inside(pOverPhi)==kOutside);
    assert(s4.Inside(pInPhi)==kInside);
    assert(s5.Inside(pbigz)==kOutside);

    assert(s41.Inside(pmx)==kSurface);
    assert(s42.Inside(pmx)==kSurface);
    assert(sn11.Inside(pzero)==kSurface);
    assert(sn12.Inside(pzero)==kSurface);
    //std::cout<<"sn11.Inside0="<<sn11.Inside(pzero)<<std::endl;
    //std::cout<<"sn12.Inside0="<<sn12.Inside(pzero)<<std::endl;
    // std::cout<<"sn13.Inside0="<<sn13.Inside(pzero)<<std::endl;
    //std::cout<<"sn13.InsideZ="<<sn13.Inside(G4ThreeVector(0,0,1.))<<std::endl;
    assert(sn13.Inside(pzero)==kSurface);
    assert(sn13.Inside(-pz)==kInside);
    assert(sn12.Inside(pz)==kInside);
    assert(sn12.Inside(px)==kOutside);
    assert(sn11.Inside(px)==kSurface);
  
    //OpticalEscape
    //std::cout<<"OpticalEscape In="<<testTheta.Inside(G4ThreeVector(3.641752094329931e-06,-2.597111833013699e-05,  -1.514111263323237e-05))<<std::endl;




// Checking G4Sphere::SurfaceNormal
    G4double p2=1./std::sqrt(2.),p3=1./std::sqrt(3.);


    norm=sn1.SurfaceNormal(G4ThreeVector(0.,0.,50.));
    assert(ApproxEqual(norm,G4ThreeVector(p3,p3,p3)));
    norm=sn1.SurfaceNormal(G4ThreeVector(0.,0.,0.));
    assert(ApproxEqual(norm,G4ThreeVector(p2,p2,0.)));

    norm=sn11.SurfaceNormal(G4ThreeVector(0.,0.,0.));
    assert(ApproxEqual(norm,G4ThreeVector(0.,0.,-1.)));
    norm=sn12.SurfaceNormal(G4ThreeVector(0.,0.,0.));
    assert(ApproxEqual(norm,G4ThreeVector(0.,0.,-1.)));
    norm=sn13.SurfaceNormal(G4ThreeVector(0.,0.,0.));
    assert(ApproxEqual(norm,G4ThreeVector(0.,0.,1.)));

    norm=sn2.SurfaceNormal(G4ThreeVector(-45.,0.,0.));
    assert(ApproxEqual(norm,G4ThreeVector(p2,-p2,0.)));



    norm=s1.SurfaceNormal(ponrmax1);
    assert(ApproxEqual(norm,vx));

// Checking G4Sphere::DistanceToOut(P)
    //Full Sphere
    Dist=s1.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,50));
    Dist=s1.DistanceToOut(ponrmax1);
    assert(ApproxEqual(Dist,0));
    Dist=s1.DistanceToOut(px);
    assert(ApproxEqual(Dist,20));
    //Full Sphere with Rmin
    Dist=s2.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,0));
    Dist=s2.DistanceToOut(ponrmax1);
    assert(ApproxEqual(Dist,0));
    Dist=s2.DistanceToOut(ponrmin1);
    assert(ApproxEqual(Dist,0));
    //Sphere with full Phi and Theta section
    Dist=sn11.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,0));
    Dist=sn11.DistanceToOut(px);
    //G4cout<<"Dist=sn11.DistanceToOut(px) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0));
    Dist=sn12.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,0));
    Dist=sn12.DistanceToOut(pz);
    //G4cout<<"Dist=sn12.Inside(pz) = "<<sn12.Inside(pz)<<G4endl;
    assert(ApproxEqual(Dist,20));
    Dist=sn13.DistanceToOut(-pz);
    assert(ApproxEqual(Dist,20));
    Dist=s3.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,0));
    //G4cout<<"Dist=sn3.pzero = "<<Dist<<G4endl;

    Dist=s9.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,0));
    Dist=s9.DistanceToOut(px);
    //G4cout<<"Dist=sn11.DistanceToOut(px) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0));
    //Sphere with Phi section
    Dist=s41.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,0));
     Dist=s41.DistanceToOut(ponrmax1);
    assert(ApproxEqual(Dist,0));


// Checking G4Sphere::DistanceToOut(p,v)


        Dist=s1.DistanceToOut(pz,vz,calcNorm,pgoodNorm,pNorm);
	// G4cout<<"Dist=s1.DistanceToOut(pz,vz) = "<<Dist<<G4endl;
        assert(ApproxEqual(Dist,20.));

     Dist=s1.DistanceToOut(ponrmax1,vx,calcNorm,pgoodNorm,pNorm);
     *pNorm=pNorm->unit();
     assert(ApproxEqual(Dist,0)&&*pgoodNorm&&ApproxEqual(*pNorm,vx));

     Dist=s2.DistanceToOut(ponrmin1,vx,calcNorm,pgoodNorm,pNorm);
     *pNorm=pNorm->unit();
     assert(ApproxEqual(Dist,5)&&*pgoodNorm&&ApproxEqual(*pNorm,vx));

     Dist=s2.DistanceToOut(ponrmax2,vx,calcNorm,pgoodNorm,pNorm);
     assert(ApproxEqual(Dist,0)&&*pgoodNorm&&ApproxEqual(*pNorm,vxy));

    Dist=s1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vx)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmx)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vy)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmy)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vz)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmz)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vxy)&&*pgoodNorm);

    Dist=s4.DistanceToOut(ponphi1,vx,calcNorm,pgoodNorm,pNorm);
    //assert(ApproxEqual(Dist,0)&&ApproxEqual(pNorm->unit(),vmxmy)&&*pgoodNorm);
    Dist=s4.DistanceToOut(ponphi2,vx,calcNorm,pgoodNorm,pNorm);
    // assert(ApproxEqual(Dist,0)&&ApproxEqual(pNorm->unit(),vmxy)&&*pgoodNorm);
    Dist=s3.DistanceToOut(pontheta1,vz,calcNorm,pgoodNorm,pNorm);
    // assert(ApproxEqual(Dist,0)&&ApproxEqual(pNorm->unit(),vy)&&*pgoodNorm);
    Dist=s32.DistanceToOut(pontheta1,vmz,calcNorm,pgoodNorm,pNorm);
    //assert(ApproxEqual(Dist,0)&&ApproxEqual(pNorm->unit(),vmy)&&*pgoodNorm);
    Dist=s32.DistanceToOut(pontheta1,vz,calcNorm,pgoodNorm,pNorm);
    //assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vz)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmz)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vxy)&&*pgoodNorm);
    
    
        Dist=s2.DistanceToOut(ponrmin1,vxy,calcNorm,pgoodNorm,pNorm);
	//        G4cout<<"Dist=s2.DistanceToOut(pormin1,vxy) = "<<Dist<<G4endl;

    Dist=s2.DistanceToOut(ponrmax1,vmx,calcNorm,pgoodNorm,pNorm);
    //    G4cout<<"Dist=s2.DistanceToOut(ponxside,vmx) = "<<Dist<<G4endl;
    Dist=s2.DistanceToOut(ponrmax1,vmxmy,calcNorm,pgoodNorm,pNorm);
    //    G4cout<<"Dist=s2.DistanceToOut(ponxside,vmxmy) = "<<Dist<<G4endl;
    Dist=s2.DistanceToOut(ponrmax1,vz,calcNorm,pgoodNorm,pNorm);
    //    G4cout<<"Dist=s2.DistanceToOut(ponxside,vz) = "<<Dist<<G4endl;

    //    Dist=s2.DistanceToOut(pbigx,vx,calcNorm,pgoodNorm,pNorm);
    //    G4cout<<"Dist=s2.DistanceToOut(pbigx,vx) = "<<Dist<<G4endl;
    //    Dist=s2.DistanceToOut(pbigx,vxy,calcNorm,pgoodNorm,pNorm);
    //    G4cout<<"Dist=s2.DistanceToOut(pbigx,vxy) = "<<Dist<<G4endl;
    //    Dist=s2.DistanceToOut(pbigx,vz,calcNorm,pgoodNorm,pNorm);
    //    G4cout<<"Dist=s2.DistanceToOut(pbigx,vz) = "<<Dist<<G4endl;
    //Test Distance for phi section
     Dist=sn22.DistanceToOut(G4ThreeVector(0.,49.,0.),vmy,calcNorm,pgoodNorm,pNorm);
     assert(ApproxEqual(Dist,49.));
     Dist=sn22.DistanceToOut(G4ThreeVector(-45.,0.,0.),vx,calcNorm,pgoodNorm,pNorm);
     assert(ApproxEqual(Dist,45.));
G4cout<<"Dist from Center ="<<sn22.DistanceToOut(G4ThreeVector(0.,49.,0),G4ThreeVector(0,-1,0))<<G4endl;
G4cout<<"Dist from Center ="<<sn22.DistanceToOut(G4ThreeVector(-45.,0.,0),G4ThreeVector(1,0,0))<<G4endl;

     //
    Dist=b216.DistanceToOut(p216,v216,calcNorm,pgoodNorm,pNorm);
    
    // G4cout<<"b216.DistanceToOut(p216,v216,... = "<<Dist<<G4endl;

    // call from outside
    //    Dist=sAlex.DistanceToOut(psAlex,vsAlex,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sAlex.DistanceToOut(psAlex,vsAlex,... = "<<Dist<<G4endl;


    // Dist=s9.DistanceToIn(s9p,s9v);
    // G4cout<<"s9.DistanceToIn(s9p,s9v,... = "<<Dist<<G4endl;
    Dist=s13.DistanceToOut(G4ThreeVector(20.,0.,0.),vz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"s13.DistanceToOut(G4ThreeVector(20.,0.,0.),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,34.641016151377549));

    Dist=s13.DistanceToOut(G4ThreeVector(20.,0.,0.),vmz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"s13.DistanceToOut(G4ThreeVector(20.,0.,0.),vmz... = "<<Dist<<G4endl;
     assert(ApproxEqual(Dist,11.547005383792508));

    Dist=s14.DistanceToOut(G4ThreeVector(20.,0.,0.),vz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"s14.DistanceToOut(G4ThreeVector(20.,0.,0.),vz... = "<<Dist<<G4endl;
     assert(ApproxEqual(Dist,11.547005383792508));

    Dist=s14.DistanceToOut(G4ThreeVector(20.,0.,0.),vmz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"s14.DistanceToOut(G4ThreeVector(20.,0.,0.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,34.641016151377549));

    Dist=sb830.DistanceToOut(pb830,vb830,calcNorm,pgoodNorm,pNorm);
    G4cout<<"sb830.DistanceToOut(pb830,vb830... = "<<Dist<<G4endl;
    inside = sb830.Inside(pb830+Dist*vb830);
    G4cout<<"sb830.Inside(pb830+Dist*vb830) = "<<OutputInside(inside)<<G4endl ;

    Dist=sn11.DistanceToOut( G4ThreeVector(0.,0.,20.), vmz, calcNorm, pgoodNorm, pNorm);
    // G4cout<<"sn11.DistanceToOut(G4ThreeVector(0.,0.,20.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,20.));

    Dist=sn11.DistanceToOut(G4ThreeVector(0.,0.,0.),vmz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn11.DistanceToOut(G4ThreeVector(0.,0.,0.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));

    Dist=sn11.DistanceToOut(G4ThreeVector(0.,0.,0.),vz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn11.DistanceToOut(G4ThreeVector(0.,0.,0.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,50.));

    Dist=sn11.DistanceToOut(G4ThreeVector(10.,0.,0.),vmz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn11.DistanceToOut(G4ThreeVector(10.,0.,0.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));
   

    Dist=sn12.DistanceToOut(G4ThreeVector(0.,0.,20.),vxmz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn12.DistanceToOut(G4ThreeVector(0.,0.,20.),vxmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,20./std::sqrt(2.)));

    Dist=sn12.DistanceToOut(G4ThreeVector(0.,0.,0.),vmz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn12.DistanceToOut(G4ThreeVector(0.,0.,0.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));

    Dist=sn12.DistanceToOut(G4ThreeVector(10.,0.,10.),vz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn12.DistanceToOut(G4ThreeVector(10.,0.,10.),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,38.989794855663561179));


    Dist=sn12.DistanceToOut(G4ThreeVector(10.,0.,10.),vmz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn12.DistanceToOut(G4ThreeVector(10.,0.,10.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));


    Dist=sn14.DistanceToOut(G4ThreeVector(10.,0.,10.),vz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn14.DistanceToOut(G4ThreeVector(10.,0.,10.),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));


    Dist=sn14.DistanceToOut(G4ThreeVector(10.,0.,5.),vz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn14.DistanceToOut(G4ThreeVector(10.,0.,5.),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,5.));


    Dist=sn14.DistanceToOut(G4ThreeVector(10.,0.,5.),vmxz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn14.DistanceToOut(G4ThreeVector(10.,0.,5.),vmxz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,3.5355339059327381968));

    Dist=sn14.DistanceToOut(G4ThreeVector(10.,0.,10.),vmxz,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn14.DistanceToOut(G4ThreeVector(10.,0.,10.),vmxz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));

     Dist=sn12.DistanceToOut(G4ThreeVector(0.,0.,0.),vx,calcNorm,pgoodNorm,pNorm);
     G4cout<<"sn12.DistanceToOut(G4ThreeVector(0.,0.,0.),vx... = "<<Dist<<G4endl;
     //    assert(ApproxEqual(Dist,0.));



    iMax = 10000;
    rMax = 49.;
    zP   = -1.;

    for( i = 0; i < iMax; i++)
    {
      rRand = rMax*G4UniformRand();
      phi   = twopi*G4UniformRand();
      pRand = G4ThreeVector( rRand*std::cos(phi), rRand*std::sin(phi), zP);

      cosTheta =  0.99*G4UniformRand();

      sinTheta = std::sqrt((1.-cosTheta)*(1.+cosTheta));     
      phi   = twopi*G4UniformRand();
      vRand = G4ThreeVector( sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);

      Dist = sn15.DistanceToOut(pRand,vRand,calcNorm,pgoodNorm,pNorm);

      if( Dist*cosTheta > -2.*zP )
      {
        G4cout<<"sn15.DistanceToOut(pRand,vRand,..., i = "<<i<<G4endl;
        G4cout<<pRand.x()<<", "<<pRand.y()<<", "<<pRand.z()<<G4endl;
        G4cout<<vRand.x()<<", "<<vRand.y()<<", "<<vRand.z()<<G4endl;
	break;
      }
    }
    pRand = G4ThreeVector( -9.3576454272980740257,  10.436745988128407703,   -1.);
    vRand = G4ThreeVector( -0.064906521070262901407, 0.35249758429187183495,  0.93355910181999202102);
    Dist = sn15.DistanceToOut(pRand,vRand,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn15.DistanceToOut(pRand,vRand... = "<<Dist<<G4endl;  
    assert(ApproxEqual(Dist,1.3409673565670394701));


    zP   = 1.;

    for( i = 0; i < iMax; i++)
    {
      rRand = rMax*G4UniformRand();
      phi   = twopi*G4UniformRand();
      pRand = G4ThreeVector( rRand*std::cos(phi), rRand*std::sin(phi), zP);

      cosTheta =  -0.99*G4UniformRand();

      sinTheta = std::sqrt((1.-cosTheta)*(1.+cosTheta));     
      phi   = twopi*G4UniformRand();
      vRand = G4ThreeVector( sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);

      Dist = sn16.DistanceToOut(pRand,vRand,calcNorm,pgoodNorm,pNorm);

      if( -Dist*cosTheta > 2.*zP )
      {
        G4cout<<"sn16.DistanceToOut(pRand,vRand,..., i = "<<i<<G4endl;
        G4cout<<pRand.x()<<", "<<pRand.y()<<", "<<pRand.z()<<G4endl;
        G4cout<<vRand.x()<<", "<<vRand.y()<<", "<<vRand.z()<<G4endl;
	break;
      }
    }



    pRand = G4ThreeVector( -18.872396210750576273, 11.573660000258405134, 1);
    vRand = G4ThreeVector( -0.85674495247509219187, -0.51247617288646407641, -0.057933225631713866632);
    Dist = sn16.DistanceToOut(pRand,vRand,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"sn16.DistanceToOut(pRand,vRand... = "<<Dist<<G4endl;  
    assert(ApproxEqual(Dist,8.9849541128705734394));
   


    zP   = -2.;

    for( i = 0; i < iMax; i++)
    {
      rRand = rMax*G4UniformRand();
      phi   = twopi*G4UniformRand();

      pRand = G4ThreeVector( rRand*std::cos(phi)*std::sin(91.*deg), 
                             rRand*std::sin(phi)*std::sin(91.*deg), 
                             rRand*std::cos(91.*deg));

      cosTheta =  std::cos( 180.*deg-89.*deg*G4UniformRand() );

      sinTheta = std::sqrt((1.-cosTheta)*(1.+cosTheta));     
      phi   = twopi*G4UniformRand();
      vRand = G4ThreeVector( sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);

      Dist = sn17.DistanceToOut(pRand,vRand,calcNorm,pgoodNorm,pNorm);

      if( -Dist*cosTheta > 50. )
      {
        G4cout<<"sn17.DistanceToOut(pRand,vRand,..., i = "<<i<<G4endl;
        G4cout<<pRand.x()<<", "<<pRand.y()<<", "<<pRand.z()<<G4endl;
        G4cout<<vRand.x()<<", "<<vRand.y()<<", "<<vRand.z()<<G4endl;
	break;
      }
    }

    pRand = G4ThreeVector(16.20075504802145, -22.42917454903122, -41.6469406430184);
    vRand = G4ThreeVector( -0.280469198715188, 0.4463870649961534, 0.8497503261406731 );
    norm = sn17.SurfaceNormal(pRand);
    G4cout<<"norm = sn17.SurfaceNormal(pRand) = "<<norm.x()<<", "<<norm.y()<<", "<<norm.z()<<G4endl;
    inside = sn17.Inside(pRand);
    G4cout<<"sn17.Inside(pRand) = "<<OutputInside(inside)<<G4endl ;
    Dist = sn17.DistanceToOut(pRand,vRand,calcNorm,pgoodNorm,pNorm);    
    G4cout<<"sn17.DistanceToOut(pRand,vRand,...) = "<<Dist<<G4endl;  
    // assert(ApproxEqual(Dist,8.9849541128705734394));
   
    pRand = G4ThreeVector(2.469342694407535, -0.5746362009323557, -0.04425421860130281);
    vRand = G4ThreeVector( -0.2513630044708909, 0.4396138160169732, -0.8622971255607673);
    norm = sn17.SurfaceNormal(pRand);
    G4cout<<"norm = sn17.SurfaceNormal(pRand) = "<<norm.x()<<", "<<norm.y()<<", "<<norm.z()<<G4endl;
    inside = sn17.Inside(pRand);
    G4cout<<"sn17.Inside(pRand) = "<<OutputInside(inside)<<G4endl ;
    Dist = sn17.DistanceToOut(pRand,vRand,calcNorm,pgoodNorm,pNorm);    
    G4cout<<"sn17.DistanceToOut(pRand,vRand,...) = "<<Dist<<G4endl;  
    Dist = sn17.DistanceToIn(pRand,vRand);    
    G4cout<<"sn17.DistanceToIn(pRand,vRand) = "<<Dist<<G4endl;  
    // assert(ApproxEqual(Dist,8.9849541128705734394));
   



// Checking G4Sphere::DistanceToIn(P)

    Dist=s2.DistanceToIn(pzero);
    assert(ApproxEqual(Dist,45));
    Dist=s1.DistanceToIn(ponrmax1);
    assert(ApproxEqual(Dist,0));

// Checking G4Sphere::DistanceToIn(P,V)


    zP   = 3.;

    for( i = 0; i < iMax; i++)
    {
      rRand = 0.5*rMax*G4UniformRand();
      phi   = twopi*G4UniformRand();
      pRand = G4ThreeVector( rRand*std::cos(phi), rRand*std::sin(phi), zP);

      cosTheta =  std::cos( 180.*deg - 78.*deg*G4UniformRand() );
      cosTheta =  -0.99; // std::cos( 180.*deg - 78.*deg*G4UniformRand() );

      sinTheta = std::sqrt((1.-cosTheta)*(1.+cosTheta));     
      phi   = twopi*G4UniformRand();
      vRand = G4ThreeVector( sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);

      Dist = sn17.DistanceToIn(pRand,vRand);

      if( -Dist*cosTheta > 2.*zP )
      {
        G4cout<<"sn17.DistanceToIn(pRand,vRand), i = "<<i<<G4endl;
        G4cout<<pRand.x()<<", "<<pRand.y()<<", "<<pRand.z()<<G4endl;
        G4cout<<vRand.x()<<", "<<vRand.y()<<", "<<vRand.z()<<G4endl;
	break;
      }
    }



    pRand = G4ThreeVector( -10.276604989981144911, -4.7072500741022338389, 1);
    vRand = G4ThreeVector( -0.28873162920537848164, 0.51647890054938394577, -0.80615357816219324061);
    Dist = sn17.DistanceToIn(pRand,vRand);
    G4cout<<"sn17.DistanceToIn(pRand,vRand) = "<<Dist<<G4endl;  
    // assert(ApproxEqual(Dist,8.9849541128705734394));
   



    Dist=s1.DistanceToIn(ponzmax,vz);
    G4cout<<"s1.DistanceToIn(ponzmax,vz) = "<<Dist<<G4endl;
    Dist=s1.DistanceToIn(pbigy,vy);
    assert(Dist==kInfinity);
    Dist=s1.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(Dist,50));

    Dist=s2.DistanceToIn(pzero,vy);
     G4cout<<"s2.DistanceToIn(pzero,vx) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,45));
    Dist=s2.DistanceToIn(pzero,vmy);
    assert(ApproxEqual(Dist,45));
    Dist=s2.DistanceToIn(ponrmin1,vx);
    //    G4cout<<"s2.DistanceToIn(ponmin1,vx) = "<<Dist<<G4endl;
    assert(Dist==0);
    Dist=s2.DistanceToIn(ponrmin1,vmx);
    assert(ApproxEqual(Dist,90));

    Dist=s2.DistanceToIn(ponrmin2,vx);
    assert(Dist==0);
    Dist=s2.DistanceToIn(ponrmin2,vmx);
    assert(ApproxEqual(Dist,90/std::sqrt(2.)));


    Dist=s3.DistanceToIn(ptesttheta1,vmz);
    //    G4cout<<"s3.DistanceToIn(ptesttheta1,vmz) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,100-48/std::sqrt(2.)));
    Dist=s3.DistanceToIn(pontheta1,vz);
    //    G4cout<<"s3.DistanceToIn(pontheta1,vz) = "<<Dist<<G4endl;
    assert(Dist==kInfinity);
    Dist=s3.DistanceToIn(pontheta1,vmz);
    assert(Dist==0);
    Dist=s3.DistanceToIn(pontheta2,vz);
    //    G4cout<<"s3.DistanceToIn(pontheta2,vz) = "<<Dist<<G4endl;
    assert(Dist==0);
    Dist=s3.DistanceToIn(pontheta2,vmz);
    assert(Dist==kInfinity);
    Dist=s32.DistanceToIn(pontheta1,vz);
    //    G4cout<<"s32.DistanceToIn(pontheta1,vz) = "<<Dist<<G4endl;
    assert(Dist==0);
    Dist=s32.DistanceToIn(pontheta1,vmz);
    //    G4cout<<"s32.DistanceToIn(pontheta1,vmz) = "<<Dist<<G4endl;
    assert(Dist==kInfinity);
    Dist=s33.DistanceToIn(pontheta2,vz);
    //    G4cout<<"s33.DistanceToIn(pontheta2,vz) = "<<Dist<<G4endl;
    assert(Dist==kInfinity);
    Dist=s33.DistanceToIn(pontheta2,vmz);
    //    G4cout<<"s33.DistanceToIn(pontheta2,vmz) = "<<Dist<<G4endl;
    assert(Dist==0);

     Dist=s4.DistanceToIn(pbigy,vmy);
     assert(Dist==kInfinity);
     Dist=s4.DistanceToIn(pbigz,vmz);
     assert(ApproxEqual(Dist,50));
     Dist=s4.DistanceToIn(pzero,vy);
     assert(Dist==kInfinity);
     Dist=s4.DistanceToIn(pzero,vx);
     assert(ApproxEqual(Dist,45));

     Dist=s4.DistanceToIn(ptestphi1,vx);
     assert(ApproxEqual(Dist,100+45/std::sqrt(2.)));
     Dist=s4.DistanceToIn(ponphi1,vmxmy);
     assert(Dist==kInfinity);
     Dist=s4.DistanceToIn(ponphi1,vxy);
     //     G4cout<<"s4.DistanceToIn(ponphi1,vxy) = "<<Dist<<G4endl;
     assert(ApproxEqual(Dist,0));

     Dist=s4.DistanceToIn(ptestphi2,vx);
     assert(ApproxEqual(Dist,100+45/std::sqrt(2.)));
     Dist=s4.DistanceToIn(ponphi2,vmxy);
     assert(Dist==kInfinity);
     Dist=s4.DistanceToIn(ponphi2,vxmy);
     //     G4cout<<"s4.DistanceToIn(ponphi2,vxmy) = "<<Dist<<G4endl;
     assert(ApproxEqual(Dist,0));

     Dist=s3.DistanceToIn(pzero,vx);
     assert(ApproxEqual(Dist,45));
     Dist=s3.DistanceToIn(ptesttheta1,vmz);
     assert(ApproxEqual(Dist,100-48/std::sqrt(2.)));
     Dist=b216.DistanceToIn(p216,v216);
     // G4cout<<"b216.DistanceToIn(p216,v216) = "<<Dist<<G4endl;

     Dist=b658.DistanceToIn(pzero,vz);
     G4cout<<"b658.DistanceToIn(pzero,vz) = "<<Dist<<G4endl;



    Dist=s13.DistanceToIn(G4ThreeVector(20.,0.,-70.),vz);
    // G4cout<<"s13.DistanceToIn(G4ThreeVector(20.,0.,-70.),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,58.452994616207498));

    Dist=s13.DistanceToIn(G4ThreeVector(20.,0.,70.),vmz);
    // G4cout<<"s13.DistanceToIn(G4ThreeVector(20.,0.,70.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,35.358983848622451));


    Dist=s14.DistanceToIn(G4ThreeVector(20.,0.,-70.),vz);
    // G4cout<<"s14.DistanceToIn(G4ThreeVector(20.,0.,-70.),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,35.358983848622451));

    Dist=s14.DistanceToIn(G4ThreeVector(20.,0.,70.),vmz);
    // G4cout<<"s14.DistanceToIn(G4ThreeVector(20.,0.,70.),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,58.452994616207498));

    Dist=b1046.DistanceToIn(G4ThreeVector(0.,0.,4800*km),vmz);
    G4cout<<"b1046.DistanceToIn(G4ThreeVector(0.,0.,4800*km),vmz... = "<<Dist<<G4endl;
    // assert(ApproxEqual(Dist,0.));
//Test for Sphere with Theta section(sn11) for point(0,0,0)
    Dist=sn11.DistanceToIn(G4ThreeVector(0.,0.,0.),vmz);
    //G4cout<<"sn11.DistanceToIn(G4ThreeVector(0.,0.,0),vmz... = "<<Dist<<G4endl; 
    assert(ApproxEqual(Dist,kInfinity));
    Dist=sn11.DistanceToIn(G4ThreeVector(0.,0.,0.),vz);
    // G4cout<<"sn11.DistanceToIn(G4ThreeVector(0.,0.,0),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));

    Dist=sn11.DistanceToOut(G4ThreeVector(0.,0.,0.),vmz);
    //G4cout<<"sn11.DistanceToOut(G4ThreeVector(0.,0.,0),vmz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.));
    Dist=sn11.DistanceToOut(G4ThreeVector(0.,0.,0.),vz);
    //G4cout<<"sn11.DistanceToIn(G4ThreeVector(0.,0.,0),vz... = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,50.));
    
    Dist=sn11.DistanceToIn(G4ThreeVector(0.0,0.,0.),vmxy);
    G4cout<<"sn11.DistanceToIn(G4ThreeVector(0.,0.,0),vmz... = "<<Dist<<G4endl; 
    //assert(ApproxEqual(Dist,kInfinity));
    Dist=sn11.DistanceToIn(G4ThreeVector(0.0,0.,0.),vxy);
    G4cout<<"sn11.DistanceToIn(G4ThreeVector(0.,0.,0),vz... = "<<Dist<<G4endl;
    //assert(ApproxEqual(Dist,0.));

     Dist=sn11.DistanceToOut(G4ThreeVector(0.0,0.,0.),vmxy);
     G4cout<<"sn11.DistanceToOut(G4ThreeVector(0.,0.,0),vmz... = "<<Dist<<G4endl;
     //assert(ApproxEqual(Dist,0.));
    Dist=sn11.DistanceToOut(G4ThreeVector(0.0,0.,0.),vxy);
    G4cout<<"sn11.DistanceToOut(G4ThreeVector(0.,0.,0),vz... = "<<Dist<<G4endl;
    //assert(ApproxEqual(Dist,50.));


    Dist=sn12.DistanceToIn(G4ThreeVector(0.0,0.,0.),vx);
    G4cout<<"sn12.DistanceToIn(G4ThreeVector(0.,0.,0),vx... = "<<Dist<<G4endl; 
    //assert(ApproxEqual(Dist,kInfinity));
    Dist=sn12.DistanceToIn(G4ThreeVector(0.0,0.,0.),vy);
    G4cout<<"sn12.DistanceToIn(G4ThreeVector(0.,0.,0),vy... = "<<Dist<<G4endl;
    //assert(ApproxEqual(Dist,0.));

     Dist=sn12.DistanceToOut(G4ThreeVector(0.0,0.,0.),vx);
     G4cout<<"sn12.DistanceToOut(G4ThreeVector(0.,0.,0),vx... = "<<Dist<<G4endl;
     //assert(ApproxEqual(Dist,0.));
    Dist=sn12.DistanceToOut(G4ThreeVector(0.0,0.,0.),vy);
    G4cout<<"sn12.DistanceToOut(G4ThreeVector(0.,0.,0),vy... = "<<Dist<<G4endl;
    //assert(ApproxEqual(Dist,50.));

    Dist=sn12.DistanceToIn(G4ThreeVector(0.0,0.,0.),vmz);
    G4cout<<"sn12.DistanceToIn(G4ThreeVector(0.,0.,0),vmz... = "<<Dist<<G4endl; 
    //assert(ApproxEqual(Dist,kInfinity));
    Dist=sn12.DistanceToIn(G4ThreeVector(0.0,0.,0.),vz);
    G4cout<<"sn12.DistanceToIn(G4ThreeVector(0.,0.,0),vz... = "<<Dist<<G4endl;
    //assert(ApproxEqual(Dist,0.));

  

     ///////////////////////////////////////////////////////////////////////////


// CalculateExtent
    G4VoxelLimits limit;		// Unlimited
    G4RotationMatrix noRot;
    G4AffineTransform origin;
    G4double min,max;
    assert(s1.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);

    assert(s2.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s2.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s2.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);

    assert(s3.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s3.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s3.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(min<=-50/std::sqrt(2.)&&max>=50/std::sqrt(2.));
    
    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);
    assert(s1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(min<=-150&&max>=-50);
    assert(s1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(min<=-160&&max>=-60);
    assert(s1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(min<=-170&&max>=-70);

    assert(s3.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(min<=-150&&max>=-50);
    assert(s3.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(min<=-160&&max>=-60);
    assert(s3.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(min<=-170+50/std::sqrt(2.)&&max>=-70-50/std::sqrt(2.));
    
    G4RotationMatrix r90Z;
    r90Z.rotateZ(halfpi);
    G4AffineTransform tRotZ(r90Z,pzero);
    assert(s1.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);

    assert(s2.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s2.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s2.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);

    assert(s3.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s3.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s3.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(min<=-50/std::sqrt(2.)&&max>=50/std::sqrt(2.));
    
// Check that clipped away
    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-60);
    assert(!s1.CalculateExtent(kXAxis,xClip,origin,min,max));

// Assert clipped to volume
    G4VoxelLimits allClip;
    allClip.AddLimit(kXAxis,-5,+5);
    allClip.AddLimit(kYAxis,-5,+5);
    allClip.AddLimit(kZAxis,-5,+5);
    G4RotationMatrix genRot;
    genRot.rotateX(pi/6);
    genRot.rotateY(pi/6);
    genRot.rotateZ(pi/6);
    G4AffineTransform tGen(genRot,vx);
    assert(s1.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s1.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s1.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);

    assert(s2.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s2.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s2.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    
    s34.CalculateExtent(kXAxis,allClip,tGen,min,max);
    // G4cout<<"s34.CalculateExtent(kXAxis,allClip,tGen,min,max)"<<G4endl ;
    // G4cout<<"min = "<<min<<"   max = "<<max<<G4endl ;

    s34.CalculateExtent(kYAxis,allClip,tGen,min,max);
    // G4cout<<"s3.CalculateExtent(kYAxis,allClip,tGen,min,max)"<<G4endl ;
    // G4cout<<"min = "<<min<<"   max = "<<max<<G4endl ;

    s34.CalculateExtent(kZAxis,allClip,tGen,min,max);
    //  G4cout<<"s34.CalculateExtent(kZAxis,allClip,tGen,min,max)"<<G4endl ;
    //  G4cout<<"min = "<<min<<"   max = "<<max<<G4endl ;
    assert(s34.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s34.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s34.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);

// Test z clipping ok
    for (G4double zTest=-100;zTest<100;zTest+=9)
	{
	    G4VoxelLimits zTestClip;
	    zTestClip.AddLimit(kZAxis,-kInfinity,zTest);
	    if (zTest<-50)
		{
		    assert(!s1.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		    assert(!s2.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		}
	    else
		{
		    assert(s1.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		    assert(s2.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		    G4double testMin=-50;
		    G4double testMax=(zTest<50) ? zTest : 50;
		    assert (ApproxEqual(min,testMin)
			    &&ApproxEqual(max,testMax));
		}
	}

// Test y clipping ok
    for (G4double xTest=-100;xTest<100;xTest+=9)
	{
	    G4VoxelLimits xTestClip;
	    xTestClip.AddLimit(kXAxis,-kInfinity,xTest);
	    if (xTest<-50)
		{
		    assert(!s1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
		}
	    else
		{
		   assert(s1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
// Calc max y coordinate
		   G4double testMax=(xTest<0) ? std::sqrt(50*50-xTest*xTest) : 50;
		   assert (ApproxEqual(min,-testMax)
			   &&ApproxEqual(max,testMax));
		}
	}

// Test x clipping ok
    for (G4double yTest=-100;yTest<100;yTest+=9)
	{
	    G4VoxelLimits yTestClip;
	    yTestClip.AddLimit(kYAxis,-kInfinity,yTest);
	    if (yTest<-50)
		{
		    assert(!s1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
		}
	    else
		{
		   assert(s1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
// Calc max y coordinate
		   G4double testMax=(yTest<0) ? std::sqrt(50*50-yTest*yTest) : 50;
		   assert (ApproxEqual(min,-testMax)
			   &&ApproxEqual(max,testMax));
		}
	}

        G4bool checkPoint( const G4Sphere& pSph, G4ThreeVector origin,
                           G4double  d,    G4ThreeVector dir,    EInside  exp); 

  G4Sphere SpAroundX("SpAroundX",  10.*mm, 1000.*mm, -1.0*degree, 2.0*degree, 0.*degree, 180.0*degree );

	G4double  sinOneDeg = std::sin( 1.0 * degree );
	  radOne = 100.0 * mm;

	G4ThreeVector  ptPhiSurfExct= G4ThreeVector( radOne * std::cos( -1.0 * degree ) , 
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
        G4ThreeVector  StartPt( radOne * std::cos(phiPoint), radOne * std::sin(phiPoint), 0.0); 
        G4cout << "For sphere " << SphDeepNeg.GetName() << G4endl;
        G4cout << " Starting from point " << ptPhiSurfExct << G4endl;

        checkPoint( SphDeepNeg, StartPt,  0.0,  vy,   kInside); 

        // Try the edges  
  G4ThreeVector  NegEdgePt( radOne * std::cos(-270.0*degree), radOne * std::sin(-270.0*degree), 0.0); 
        G4ThreeVector  PosEdgePt( radOne * std::cos(10.0*degree), radOne * std::sin(10.0*degree), 0.0); 

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



	// G4double sTheta = 135*degree, pxt = -10., pyt = 10.;

	// G4cout<<"sTheta = "<<sTheta/degree<<" degree"<<G4endl;
        // G4cout<<"std::tan(sTheta) = "<<std::tan(sTheta)<<G4endl;
        // G4cout<<"std::sin(sTheta) = "<<std::sin(sTheta)<<G4endl;
        // G4cout<<"std::atan(sTheta) = "<<std::atan(sTheta)/degree<<G4endl;
        // G4cout<<"std::atan2(pyt,pxt) = pyt/pxt = "<<std::atan2(pyt,pxt)/degree<<G4endl;

        //Check center point
   

         G4Sphere sntest("sntest",0,50,0,twopi,0.0,0.75*pi);
         G4ThreeVector in1, in2;
         in1=G4ThreeVector(0,0,30);
         in2=G4ThreeVector(0,0,0);
         G4double length = (in1-in2).mag();
         G4ThreeVector dir12 = (in2-in1)/length;
         for(i = 0 ; i < 100; i++)
	   {
             if(sntest.Inside(in1) != kOutside)
	       {Dist=sntest.DistanceToOut(in1,dir12);
		 if(std::fabs(Dist-length)>0.00001) G4cout<<"Potential ERROR i="<<i<<" Dout="<<Dist<<" dif="<< Dist-length<<G4endl;
                in1=in1+G4ThreeVector(0.00001,0.00001,0);
                //in2=in2+G4ThreeVector(-0.00001,-0.00001,0);
                length = (in1-in2).mag();
                dir12 = (in2-in1)/length;
	       }
             else{ G4cout<<" error in test In="<<sntest.Inside(in1)<<G4endl;}
	   }


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
    G4int verbose = 0; 

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
       G4cout << "  Origin=    " << origin << G4endl; 
       G4cout << "  Direction= " << direction  << G4endl; 
       G4cout << "  dist= " << dist;
       G4cout << "  Actual-point= " << newPoint << G4endl;    
    } 

    distOut = rSphere.DistanceToOut( newPoint, direction ); 
                    /*=============*/
    if ( verbose ) G4cout << " DistToOut (p, dir) = " << distOut << G4endl;
 
    G4cout.precision(oldPrecision); 

    return good;
}

