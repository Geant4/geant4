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


#include <assert.h>
#include <cmath>

#include "globals.hh"
#include "geomdefs.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Orb.hh"
#include "G4Ellipsoid.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"

#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"


// #include "G4DisplacedSolid.hh"

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


int main()
{
  G4int i;
  G4ThreeVector pzero(0,0,0), p;

  G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40),

                    ponb2x(10,0,0),  ponb2y(0,10,0),  ponb2z(0,0,10),

                   ponb2mx(-10,0,0),ponb2my(0,-10,0),ponb2mz(0,0,-10);

  G4ThreeVector ponmxside(-20,0,0),ponmyside(0,-30,0),ponmzside(0,0,-40);

  G4ThreeVector ponzsidey(0,25,40),ponmzsidey(0,25,-40),

                  ponb2zy(0,5,10),ponb2mzy(0,5,-10) ;

  G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);

  G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

  G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);

  G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);

  G4ThreeVector vxy(1/std::sqrt(2.0),1/std::sqrt(2.0),0);

  G4ThreeVector vmxy(-1/std::sqrt(2.0),1/std::sqrt(2.0),0);

  G4ThreeVector vmxmy(-1/std::sqrt(2.0),-1/std::sqrt(2.0),0);

  G4ThreeVector vxmy(1/std::sqrt(2.0),-1/std::sqrt(2.0),0);

  G4ThreeVector vxmz(1/std::sqrt(2.0),0,-1/std::sqrt(2.0));

  G4double dist;
  G4ThreeVector *pNorm,norm;
  G4bool *pgoodNorm,goodNorm,calcNorm=true;

  pNorm=&norm;
  pgoodNorm=&goodNorm;

  G4RotationMatrix identity, xRot ;
    
// NOTE: xRot = rotation such that x axis->-x axis & y axis->-y axis

  xRot.rotateZ(-pi);

  G4Transform3D transform(xRot,G4ThreeVector(0,30,0)) ;
  G4Transform3D transform2(xRot,G4ThreeVector(50,0,0)) ;

  G4Box b1("Test Box #1",20.,30.,40.);
  G4Box b2("Test Box #2",10.,10.,10.);
  G4Box b3("Test Box #3",10.,20.,50.);
  G4Box b4("Test Box #4",20.,20.,40.);
  G4Box b5("Test Box #4",50.,50.,50.);

  G4Tubs t1("Solid Tube #1",0,50.,50.,0.,360.*degree);
    
    // t2\t3 for DistanceToIn

  G4Tubs t2("Hole Tube #2",50.,60.,50.,0.,2.*pi); 
 
  G4Tubs t3("Hole Tube #3",45.,55.,50.,pi/4.,pi*3./2.);

  G4Cons c1("Hollow Full Tube",50.,100.,50.,100.,50.,0,2*pi),

	   c2("Full Cone",0,50.,0,100.,50.,0,2*pi) ;

  G4SubtractionSolid b1Sb2("b1Sb2",&b1,&b2),

                       t1Sb2("t1Sb2",&t1,&b2),

                       t2St3("t2St3",&t2,&t3),

                       c2Sb2("c2Sb2",&c2,&b2) ;

    // Tube t1 \ box t3, which was rotated by -pi and translated +30y

  G4SubtractionSolid   t1Sb3("t1Subtractionb3",&t1,&b3,transform) ;
  G4SubtractionSolid   b1Sb4("t1Subtractionb3",&b1,&b4,transform) ;
  G4SubtractionSolid   b5St1("b4St1",&b5,&t1,transform2) ;

  G4SubtractionSolid   b1Sb2touch("b1Sb4touch",&b1,&b2,&identity,G4ThreeVector(10.,0,0)) ;





  G4Tubs* tube4 = new G4Tubs( "OuterFrame",
                      1.0*m,            // inner radius
                      1.1*m,            // outer radius
                      0.01*m,           // half-thickness in z
                      -15*deg,          // start angle
                      30*deg );         // total angle

  G4Box* tube5 = new G4Box( "Cutout",    // name
                   0.02*m,      // half-width (x)
                   0.25*m,      // half-height (y)
                   0.01001*m ); // half-thickness (z)

  G4Transform3D tran2 = G4Translate3D( 1.03*m, 0.0, 0.0 );

  G4VSolid* solid = new G4SubtractionSolid( "drcExample",tube4,tube5, tran2 );

  G4Cons* cone3 = new G4Cons( "OuterFrame",
                              0.6*m, // pRmin1
                              1.0*m, // pRmax1
                              0.2*m, // pRmin2
                              0.8*m, // pRmax2
                              0.2*m,
                              0*deg,
                              180*deg );

  G4Cons* cone4 = new G4Cons( "OuterFrame",
                              0.6*m, // pRmin1
                              1.0*m, // pRmax1
                              0.2*m, // pRmin2
                              0.8*m, // pRmax2
                              0.2*m,
                              0*deg,
                              180*deg );
  G4RotationMatrix rotmat3;
  rotmat3.rotateY(pi/4.0);
  G4Transform3D tran3 = G4Transform3D(rotmat3,G4ThreeVector(0.0,0.0,0.0));

  G4VSolid* c3Ic4 = new G4SubtractionSolid( "Example", cone3, cone4, tran3 );


  //G4Torus* insp    = new G4Torus("Isp", 7.5*cm, 8.1*cm, 15.6*cm,  0, 2*pi);
  //G4Tubs*  tmptube = new G4Tubs("TMPT", 7.5*cm, 15.6*cm, 8.1*cm, 0, 2*pi);

  G4RotationMatrix* rotDz5 = new G4RotationMatrix();
  rotDz5->rotateZ(pi/2.);

  // G4SubtractionSolid* inn = new G4SubtractionSolid("Innerpart", insp,
  //                                  tmptube, rotDz5, G4ThreeVector(0,0,0));

  // LHCB RICH-detector problem:

  G4ThreeVector lhcbTransl(8606.,173396.33,100110.42);
  G4RotationMatrix lhcbMat;
  lhcbMat.rotateX(-60*degree);
  G4Transform3D lhcbTran = G4Transform3D(lhcbMat,lhcbTransl);
  G4Sphere* lhcbSphere = new G4Sphere("lhcbSphere",8600*mm, 
				                   8606*mm, 
				             -1.6991355*degree, 
				              3.3982711*degree, 
				              88.528559*degree, 
				              2.9428812*degree   );

  G4Box* lhcbBox = new G4Box("lhcbBox",100000*mm, 
			               100000*mm, 
                                       200000*mm   );

  
  G4VSolid* lhcbSub = new G4SubtractionSolid( "lhcbSub", lhcbSphere, 
                                                       lhcbBox, lhcbTran );

  G4ThreeVector lhcbP(7298.73956059975,-394.9290932077643,1497.43923688271);
  G4ThreeVector lhcbV(0.6062807093945133,0.06537967965081047,-0.7925586406726276);


  G4Orb* largeOrb = new G4Orb("largeOrb",11.*mm);
  G4Orb* smallOrb = new G4Orb("smallOrb",10.*mm);

  G4ThreeVector zeroVec(0.0, 0.0, 0.0); 
  G4VSolid* orbShell = new G4SubtractionSolid( "orbShell", largeOrb, 
                                                       smallOrb, 0, zeroVec );

  G4Ellipsoid* largeEll = new G4Ellipsoid("largeEll",11.*mm, 11*mm, 11*mm); //, -11*mm, 11*mm);
  G4Ellipsoid* smallEll = new G4Ellipsoid("smallEll",10.*mm, 10*mm, 10*mm); //, -10*mm, 10*mm);

  G4VSolid* ellShell = new G4SubtractionSolid( "orbShell", largeEll, 
                                                       smallEll, 0, zeroVec );

  G4Box* orbBox = new G4Box("orbBox",20*mm, 20*mm, 20*mm   );

  G4RotationMatrix orbMat;
  G4ThreeVector orbTransl(0.,0.,-20.);
  G4Transform3D orbTran = G4Transform3D(orbMat, orbTransl);

  G4VSolid* orbMirror = new G4IntersectionSolid( "orbMirror", orbShell, 
                                                       orbBox, orbTran);

  G4VSolid* ellMirror = new G4IntersectionSolid( "orbMirror", ellShell, 
                                                       orbBox, orbTran);

   G4cout.precision(16) ;

// Check Inside

    assert(b1Sb2.Inside(pzero)==kOutside);
    assert(b1Sb2.Inside(pbigz)==kOutside);
    assert(b1Sb2.Inside(ponb2x)==kSurface);
    assert(b1Sb2.Inside(ponb2y)==kSurface);
    assert(b1Sb2.Inside(ponb2z)==kSurface);

    assert(t1Sb3.Inside(pzero)==kInside);
    assert(t1Sb3.Inside(pbigz)==kOutside);
    assert(t1Sb3.Inside(ponb2x)==kInside);
    assert(t1Sb3.Inside(ponb2y)==kSurface);
    assert(t1Sb3.Inside(ponb2z)==kInside);

    assert(c2Sb2.Inside(pzero)==kOutside);
    assert(c2Sb2.Inside(pbigz)==kOutside);
    assert(c2Sb2.Inside(ponb2x)==kSurface);
    assert(c2Sb2.Inside(ponb2y)==kSurface);
    assert(c2Sb2.Inside(ponb2z)==kSurface);

    assert(b1Sb2touch.Inside(G4ThreeVector(20.,0,0))==kOutside);
    assert(b1Sb2touch.Inside(G4ThreeVector(20.,10.,0))==kSurface);

  

// Check Surface Normal

    G4ThreeVector normal;

    normal=b1Sb2.SurfaceNormal(ponb2x);
    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));

    normal=b1Sb2.SurfaceNormal(ponb2mx);
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));

    normal=b1Sb2.SurfaceNormal(ponb2y);
    assert(ApproxEqual(normal,G4ThreeVector(0,-1,0)));

    normal=b1Sb2.SurfaceNormal(ponb2my);
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));

    normal=b1Sb2.SurfaceNormal(ponb2z);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));

    normal=b1Sb2.SurfaceNormal(ponb2mz);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));

    normal=b1Sb2.SurfaceNormal(ponb2zy);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));

    normal=b1Sb2.SurfaceNormal(ponb2mzy);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));


    normal=t1Sb3.SurfaceNormal(ponb2y);
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));


// DistanceToOut(P)

    dist=b1Sb2.DistanceToOut(ponb2x);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToOut(ponb2y);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToOut(ponb2z);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToOut(ponxside);
    assert(ApproxEqual(dist,0));

    dist=t1Sb3.DistanceToOut(ponb2y);
    assert(ApproxEqual(dist,0));

    dist=t1Sb3.DistanceToOut(pzero);
    assert(ApproxEqual(dist,10));

    dist=t1Sb3.DistanceToOut(ponb2x);
    assert(ApproxEqual(dist,10));

    dist=t1Sb3.DistanceToOut(ponb2z);
    assert(ApproxEqual(dist,10));


// DistanceToOut(P,V)

    dist=b1Sb2.DistanceToOut(ponb2x,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponxside,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmx)); // &&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2y,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponyside,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vmy)); // &&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2z,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vz)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponzside,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vmz)); // &&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2x,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,std::sqrt(200.))&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2x,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2x,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vmx)); // &&*pgoodNorm);

    dist=b1.DistanceToOut(ponxside,vy,calcNorm,pgoodNorm,pNorm);

//  G4cout<<"b1.DistanceToOut(ponxside,vy) = "<<dist<<G4endl;
//  assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2mx,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmx)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2y,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2my,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vmy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2z,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vz)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2mz,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vmz)&&*pgoodNorm);

    // With placement

    dist=t1Sb3.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vy)); // &&*pgoodNorm);

    dist=t1Sb3.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vx)&&*pgoodNorm);

    dist=t1Sb3.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vmy)&&*pgoodNorm);

    dist=t1Sb3.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vz)&&*pgoodNorm);


//DistanceToIn(P)

    dist=b1Sb2.DistanceToIn(pbigx);
    assert(ApproxEqual(dist,80));

    dist=b1Sb2.DistanceToIn(ponb2x);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToIn(ponxside);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToIn(pbigmy);
    assert(ApproxEqual(dist,70));

    dist=b1Sb2.DistanceToIn(pbigz);
    assert(ApproxEqual(dist,60));

    dist=b1Sb2.DistanceToIn(pbigmz);
    assert(ApproxEqual(dist,60));

    // With placement

    dist=t1Sb3.DistanceToIn(pbigy);
    assert(ApproxEqual(dist,50));

    dist=t1Sb3.DistanceToIn(pbigmy);
    assert(ApproxEqual(dist,50));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,20,0));
    assert(ApproxEqual(dist,10));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,15,0));
    assert(ApproxEqual(dist,5));

// DistanceToIn(P,V)

    dist=b1Sb2.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(dist,80));

    dist=b1Sb2.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(dist,80));

    dist=b1Sb2.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(dist,70));

    dist=b1Sb2.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(dist,70));

    dist=b1Sb2.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(dist,60));

    dist=b1Sb2.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(dist,60));

    dist=b1Sb2.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    dist=b1Sb2.DistanceToIn(pbigmx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    // Cases with a number of iterations inside 'while'

    dist=t2St3.DistanceToIn(G4ThreeVector(2.5,-52.5,0),vy) ;
    assert(ApproxEqual(dist,107.443));
// G4cout<<"t2St3.DistanceToIn(G4ThreeVector(2.5,-52.5,0),vy) = "<<dist<<G4endl ;

    dist=t2St3.DistanceToIn(G4ThreeVector(2.5,-62.5,0),vy) ;
    assert(ApproxEqual(dist,2.55211));
//    G4cout<<"t2St3.DistanceToIn(G4ThreeVector(2.5,-62.5,0),vy) = "<<dist<<G4endl ;

    // With placement

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,100,0),vmy);
    G4cout<<"t1Sb3.DistanceToIn(G4ThreeVector(0,100,0),vmy) = "<<dist<<G4endl; // 90?
    // assert(ApproxEqual(dist,90)); // was 50

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,36,0),vmy);
    assert(ApproxEqual(dist,26));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,36,0),vy);
    assert(ApproxEqual(dist,kInfinity));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,36,0),vx);
    assert(ApproxEqual(dist,10));

    dist=solid->DistanceToIn(
    G4ThreeVector(1462.639634076807,-1183.908767887379,1926.359643354072),
    G4ThreeVector(-0.1985206962900923,0.5639538694097884,-0.8015894000810041));
    //  G4cout<<"solid.DistanceToIn(p,v) = "<<dist<<G4endl ;
    assert(ApproxEqual(dist,2390.7));

    dist=solid->DistanceToIn(
    G4ThreeVector(1572.774744032383,-1704.268426706159,1922.635294938558),
    G4ThreeVector(-0.2113262591952278,0.6627445069650045,-0.7184086098191366));
    // G4cout<<"solid.DistanceToIn(p,v) = "<<dist<<G4endl ;
    assert(ApproxEqual(dist,2663.0611));

    dist=c3Ic4->DistanceToIn(
    G4ThreeVector(967.602560486646,-1211.626221342084,610.1637191085789),
    G4ThreeVector(-0.1174412716483462,0.8263033300403027,-0.5508451274885945));
    // assert(ApproxEqual(dist,kInfinity));
    G4cout<<"c3Ic4->DistanceToIn = "<<dist<<G4endl ;

    dist=lhcbSub->DistanceToIn(lhcbP,lhcbV);
    //  G4cout<<"lhcbSub->DistanceToIn(lhcbP,lhcbV) = "<<dist<<G4endl ;

    G4ThreeVector orbP = G4ThreeVector(0.,0.,-9.999);


    G4double orbTheta = 30*degree;
    G4double orbPhi   = 30*degree;

    G4ThreeVector orbV = G4ThreeVector(std::sin(orbTheta)*std::cos(orbPhi),
                                       std::sin(orbTheta)*std::sin(orbPhi),
                                       std::cos(orbTheta));


    EInside insideP = smallOrb->Inside(orbP);
    G4cout<<"smallOrb->Inside(orbP) = "<<OutputInside(insideP)<<G4endl ;
    dist = smallOrb->DistanceToOut(orbP,orbV);
    G4cout<<"smallOrb->DistanceToOut(orbP,orbV) = "<<dist<<G4endl ;

    insideP = orbShell->Inside(orbP);
    G4cout<<"orbShell->Inside(orbP) = "<<OutputInside(insideP)<<G4endl ;
    dist = orbShell->DistanceToIn(orbP,orbV);
    G4cout<<"orbShell->DistanceToIn(orbP,orbV) = "<<dist<<G4endl ;

    insideP = orbMirror->Inside(orbP);
    G4cout<<"orbMirror->Inside(orbP) = "<<OutputInside(insideP)<<G4endl ;
    dist = orbMirror->DistanceToIn(orbP,orbV);
    G4cout<<"orbMirror->DistanceToIn(orbP,orbV) = "<<dist<<G4endl ;

    insideP = ellShell->Inside(orbP);
    G4cout<<"ellShell->Inside(orbP) = "<<OutputInside(insideP)<<G4endl ;
    dist = ellShell->DistanceToIn(orbP,orbV);
    G4cout<<"ellShell->DistanceToIn(orbP,orbV) = "<<dist<<G4endl ;

    insideP = ellMirror->Inside(orbP);
    G4cout<<"ellMirror->Inside(orbP) = "<<OutputInside(insideP)<<G4endl ;
    dist = ellMirror->DistanceToIn(orbP,orbV);
    G4cout<<"ellMirror->DistanceToIn(orbP,orbV) = "<<dist<<G4endl ;

  for(i=0;i<3000;i++)
  {
    orbTheta = 45*degree*G4UniformRand();
    orbPhi   = 360*degree*G4UniformRand();
    orbV = G4ThreeVector(std::sin(orbTheta)*std::cos(orbPhi),
                                       std::sin(orbTheta)*std::sin(orbPhi),
			 std::cos(orbTheta));
    dist = ellMirror->DistanceToIn(orbP,orbV);

    if(dist != kInfinity)
    {
      G4cout<<i<<"\t"<<"ellMirror->DistanceToIn(orbP,orbV) = "<<dist<<G4endl ;
    }
  }
  G4cout<<G4endl;

    G4cout<<"Tracking functions are OK"<<G4endl ;


// Point on surface

  G4cout<<G4endl;
  //  G4cout<<"Point on surface of t1Sb3"<<G4endl<<G4endl;
  G4cout<<"Point on surface of b5St1"<<G4endl<<G4endl;
  for(i=0;i<30;i++)
  {
    // p = t1Sb3.GetPointOnSurface();
      p = b5St1.GetPointOnSurface();
      G4cout<<p.x()<<"\t\t"<<p.y()<<"\t\t"<<p.z()<<G4endl;
  }
  G4cout<<G4endl;



// CalculateExtent

    G4VoxelLimits limit;		// Unlimited
    G4RotationMatrix noRot;
    G4AffineTransform origin;
    G4double min,max;

    assert(b1.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));

    assert(b1.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(b1.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    assert(b1Sb4.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));

    assert(b1Sb4.CalculateExtent(kYAxis,limit,origin,min,max));
    //  G4cout<<"min of b1Sb4.CalculateExtent(kYAxis,limit,origin,min,max) = "
    //      <<min<<G4endl ;
    // G4cout<<"max of b1Sb4.CalculateExtent(kYAxis,limit,origin,min,max) = "
    //      <<max<<G4endl ;
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(b1Sb4.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));



    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);

    assert(b1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-120)&&ApproxEqual(max,-80));

    assert(b1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-140)&&ApproxEqual(max,-80));

    assert(b1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-160)&&ApproxEqual(max,-80));

    G4RotationMatrix r90Z;
    r90Z.rotateZ(pi/2);
    G4AffineTransform tRotZ(r90Z,pzero);

    assert(b1.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(b1.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));

    assert(b1.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

// Check that clipped away

    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-50);

    assert(!b1.CalculateExtent(kXAxis,xClip,origin,min,max));

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

    assert(b1.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    G4VoxelLimits buggyClip2;
    buggyClip2.AddLimit(kXAxis,5,15);

    assert(b1.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));

    assert(b1.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(b1.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    buggyClip2.AddLimit(kYAxis,5,15);

    assert(b1.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));

    assert(b1.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));

    assert(b1.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    G4VoxelLimits buggyClip1;
    buggyClip1.AddLimit(kXAxis,-5,+5);

    assert(b1.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(b1.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    buggyClip1.AddLimit(kYAxis,-5,+5);

    assert(b1.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    G4cout<<"CalculateExtent is OK "<<G4endl ;

  delete rotDz5;
  G4int out =0 ;
  return out ;
}








