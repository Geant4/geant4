// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// BREP solid test, create by L. Broglia, 20/10/98
// modification of old G4Gerep test
//



#include <stdio.h>
#include <math.h>
#include <fstream.h>
#include <stdlib.h>
#include "G4ios.hh" 
#include "G4Axis2Placement3D.hh"
#include "G4BREPSolid.hh"
#include "G4BREPSolidBox.hh"
#include "G4BREPSolidCylinder.hh"
#include "G4BREPSolidCone.hh"
#include "G4BREPSolidTorus.hh"
#include "G4BREPSolidSphere.hh"
#include "G4BREPSolidPCone.hh"
#include "G4BREPSolidPolyhedra.hh"
#include "G4Timer.hh"
#include "G4NISTStepReader.hh"
#include "G4BezierSurface.hh"


int main(int argc, char **argv)
{
  G4Timer timer;

  G4ThreeVector S1(0,0,0);
  G4ThreeVector S2(0,50,0);
  G4ThreeVector S3(50,0,500);
  G4ThreeVector S4(-50,0,-500);  

  G4ThreeVector D1(0,0,1);
  G4ThreeVector D2(0,0,1);
  G4ThreeVector D3(0,0,-1);
  G4ThreeVector D4(0,0,1);      

  G4ThreeVector S1b(0,0,10);
  G4ThreeVector S2b(0,50,10);
  G4ThreeVector S3b(50,0,10);
  G4ThreeVector S4b(50,0,10);  

  G4ThreeVector D1b(0,0,1);
  G4ThreeVector D2b(0,0,1);
  G4ThreeVector D3b(0,0,-1);
  G4ThreeVector D4b(0,1,0);      
 
  G4ThreeVector PCon1(150,0,5);
  G4ThreeVector PConD1(-1,0,0);
  G4ThreeVector PConD2(-1,0,0.01);
  G4ThreeVector PConD3(-1,0,-0.01);  

  G4ThreeVector tStart(19000,0,10.1);
  G4ThreeVector tDir(-1,0,0);
  G4ThreeVector tStart2(0,0,10);
  G4ThreeVector pt(0,0,50);
  G4ThreeVector pt2(100000,0,50);
  double Dist;

  double RMINVec[4];
  RMINVec[0] = 10;
  RMINVec[1] = 10;
  RMINVec[2] = 10;
  RMINVec[3] = 10;  

  double RMAXVec[4];
  RMAXVec[0] = 110;
  RMAXVec[1] = 100;
  RMAXVec[2] = 100;
  RMAXVec[3] = 110;  

  double Z_Values[4];
  Z_Values[0] = 0;
  Z_Values[1] = 10;
  Z_Values[2] = 20;
  Z_Values[3] = 30;
  


  G4cout << "\n=======     Half PGon test      ========";
  
  G4BREPSolidPolyhedra *MyPGone2 = new G4BREPSolidPolyhedra ("MyPolyhedra",
							     0            ,
							     pi           ,
							     2            ,
							     3            ,
							     0            ,
							     Z_Values     ,
							     RMINVec      ,
							     RMAXVec       );
  G4cout << "\n\nHalf Polygon created ! ";
  Dist = MyPGone2->DistanceToIn(S1, D1);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPGone2->DistanceToIn(S2, D2);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPGone2->DistanceToIn(S3, D3);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPGone2->DistanceToIn(S4, D4);
  G4cout << "\nDist to in : " << Dist ;

  Dist = MyPGone2->DistanceToOut(S1b, D1b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPGone2->DistanceToOut(S2b, D2b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPGone2->DistanceToOut(S3b, D3b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPGone2->DistanceToOut(S4b, D4b);
  G4cout << "\nDist to out : " << Dist ;



  G4cout << "\n=======     PGon test      ========";

  G4BREPSolidPolyhedra *MyPGone = new G4BREPSolidPolyhedra ("MyPolyhedra",
							    0            ,
							    2*pi         ,
							    4            ,
							    3            ,
							    0            ,
							    Z_Values     ,
							    RMINVec      ,
							    RMAXVec       );

  G4cout << "\n\nPolygon created ! ";
  Dist = MyPGone->DistanceToIn(S1, D1);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPGone->DistanceToIn(S2, D2);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPGone->DistanceToIn(S3, D3);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPGone->DistanceToIn(S4, D4);
  G4cout << "\nDist to in : " << Dist ;

  Dist = MyPGone->DistanceToOut(S1b, D1b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPGone->DistanceToOut(S2b, D2b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPGone->DistanceToOut(S3b, D3b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPGone->DistanceToOut(S4b, D4b);
  G4cout << "\nDist to out : " << Dist ;


  
  G4cout << "\n=======     PCon test      ========";

  G4BREPSolidPCone *MyPCone = new G4BREPSolidPCone ("MyPCone",
						    0        ,
						    2*pi     ,
						    3        ,
						    0        ,
						    Z_Values ,
						    RMINVec  ,
						    RMAXVec   );

  G4cout << "\n\nPCone created ! ";
  Dist = MyPCone->DistanceToIn(S1, D1);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPCone->DistanceToIn(S2, D2);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPCone->DistanceToIn(S3, D3);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPCone->DistanceToIn(S4, D4);
  G4cout << "\nDist to in : " << Dist ;

  Dist = MyPCone->DistanceToOut(S1b, D1b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPCone->DistanceToOut(S2b, D2b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPCone->DistanceToOut(S3b, D3b);
  G4cout << "\nDist to out : " << Dist ;
  Dist = MyPCone->DistanceToOut(S4b, D4b);
  G4cout << "\nDist to out : " << Dist ;      
  
  Dist = MyPCone->DistanceToIn(PCon1, PConD1);
  G4cout << "\nDist to in : " << Dist ;
  Dist = MyPCone->DistanceToIn(PCon1, PConD2);
  G4cout << "\nDist to in : " << Dist ;

  Dist = MyPCone->DistanceToIn(PCon1, PConD3);
  G4cout << "\nDist to in : " << Dist ;



  G4cout << "\n ============   Cone test   ================";

  G4BREPSolidCone *MyConBox = new G4BREPSolidCone ("MyConBox"          ,
						   G4ThreeVector(0,0,0),
						   G4ThreeVector(0,0,1),
						   G4ThreeVector(1,0,0),     
						   1000.0              ,
						   0.0                 ,
						   101.0                );

  G4cout << "\n\nCone created ! ";
  G4cout << "\nDir =  -1,0,0";
  G4cout << "\nStart 19000,1,0";
  Dist = MyConBox->DistanceToIn(tStart, tDir);
  G4cout << "\nDist to in : " << Dist;
  MyConBox->Reset();

  Dist = MyConBox->DistanceToOut(tStart2, tDir);  
  G4cout << "\nStart 0,0,0";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MyConBox->DistanceToOut(pt);  
  G4cout << "\nPoint 0,0,50";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MyConBox->DistanceToIn(pt2);  
  G4cout << "\nPoint 100000,0,50";
  G4cout << "\nDist to in: " << Dist ;



  G4cout << "\n ============   Sphere test   ================";

  G4BREPSolidSphere *MySphere = new G4BREPSolidSphere ("MySphere"          ,
						       G4ThreeVector(0,0,0),
						       G4ThreeVector(1,0,0),
						       G4ThreeVector(0,0,1),
						       100                  );

  G4cout << "\n\nSphere created ! ";
  G4cout << "\nDir =  -1,0,0";
  G4cout << "\nStart 19000,1,0";
  Dist = MySphere->DistanceToIn(tStart, tDir);
  G4cout << "\nDist to in : " << Dist ;

  Dist = MySphere->DistanceToOut(tStart2, tDir);  
  G4cout << "\nStart 0,0,0";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MySphere->DistanceToOut(pt);  
  G4cout << "\nPoint 0,0,50";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MySphere->DistanceToIn(pt2);  
  G4cout << "\nPoint 100000,0,50";
  G4cout << "\nDist to in: " << Dist << "\n";



  G4cout << "\n ============   Torus test   ================";

  G4BREPSolidTorus *MyTorBox = new G4BREPSolidTorus ("MyTorBox"          ,
						     G4ThreeVector(0,0,0),
						     G4ThreeVector(0,0,1),
						     G4ThreeVector(1,0,0),     
						     0.0                 ,
						     100.0                );

  
  G4cout << "\n\nTorus created ! ";
  G4cout << "\nDir =  -1,0,0";
  G4cout << "\nStart 19000,1,0";
  Dist = MyTorBox->DistanceToIn(tStart, tDir);
  G4cout << "\nDist to in : " << Dist;
  MyTorBox->Reset();

  Dist = MyTorBox->DistanceToOut(tStart2, tDir);  
  G4cout << "\nStart 0,0,0";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MyTorBox->DistanceToOut(pt);  
  G4cout << "\nPoint 0,0,50";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MyTorBox->DistanceToIn(pt2);  
  G4cout << "\nPoint 100000,0,50";
  G4cout << "\nDist to in: " << Dist ;



  G4cout << "\n ============   Box test   ================"; 

  G4BREPSolidBox *myCalBox = new G4BREPSolidBox 
    ( "MyBox",
      G4Point3D(-1500*cm, -1500*cm, -1000*cm),
      G4Point3D(-1500*cm, -1500*cm,  1000*cm),       
      G4Point3D(-1500*cm,  1500*cm,  1000*cm),
      G4Point3D(-1500*cm,  1500*cm, -1000*cm),
      G4Point3D( 1500*cm, -1500*cm, -1000*cm),       
      G4Point3D( 1500*cm, -1500*cm,  1000*cm),
      G4Point3D( 1500*cm,  1500*cm,  1000*cm),
      G4Point3D( 1500*cm,  1500*cm, -1000*cm) );

  G4cout << "\n\nBox created ! ";
  G4cout << "\nDir =  0,-1,0";
  G4cout << "\nStart 19000,1,0";
  Dist = myCalBox->DistanceToIn(tStart, tDir);
  G4cout << "\nDist to in : " << Dist;
  myCalBox->Reset();

  Dist = myCalBox->DistanceToOut(tStart2, tDir);  
  G4cout << "\nStart 0,0,0";
  G4cout << "\nDist to out: " << Dist;

  Dist = myCalBox->DistanceToOut(pt);  
  G4cout << "\nPoint 0,0,50";
  G4cout << "\nDist to out: " << Dist;

  Dist = myCalBox->DistanceToIn(pt2);  
  G4cout << "\nPoint 100000,0,0";
  G4cout << "\nDist to in: " << Dist;



  G4cout << "\n ============   Cylinder test   ================"; 

  G4BREPSolidCylinder *myCylBox= new G4BREPSolidCylinder("MyCylBox",
							 G4ThreeVector(0,0,0),
							 G4ThreeVector(0,0,1),
							 G4ThreeVector(1,0,0),
							 100.0               ,
							 1000.0              );

  G4cout << "\n\nCylinder created ! ";
  G4cout << "\nDir =  0,-1,0";
  G4cout << "\nStart 19000,1,0";
  Dist = myCylBox->DistanceToIn(tStart, tDir);
  G4cout << "\nDist to in : " << Dist;
  myCylBox->Reset();

  Dist = myCylBox->DistanceToOut(tStart2, tDir);  
  G4cout << "\nStart 0,0,0";
  G4cout << "\nDist to out: " << Dist;

  Dist = myCylBox->DistanceToOut(pt);  
  G4cout << "\nPoint 0,0,50";
  G4cout << "\nDist to out: " << Dist;

  Dist = myCylBox->DistanceToIn(pt2);  
  G4cout << "\nPoint 100000,0,50";
  G4cout << "\nDist to in: " << Dist;

  G4cout << endl << endl;
}

