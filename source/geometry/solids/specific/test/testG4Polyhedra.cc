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
// testG4Polyhedra.cc is derived from testG4Polycone.cc  13/10/2010
//
// Cvs version: $ Id $
// Cvs tag :  $ Name $

#include "G4Timer.hh"
#undef NDEBUG
#include <assert.h>
#include <cmath>
#include <fstream>
#include <stdlib.h>
#include "G4ios.hh" 
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Polyhedra.hh"

#include "globals.hh"
#include "G4Vector3D.hh"
#include "G4Point3D.hh"
#include "G4GeometryTolerance.hh"
G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
int main(int, char **)
{

 double RMINVec[8];
  RMINVec[0] = 30;
  RMINVec[1] = 30;
  RMINVec[2] =  0;
  RMINVec[3] =  0;
  RMINVec[4] =  0;  
  RMINVec[5] =  0;
  RMINVec[6] = 40;
  RMINVec[7] = 40;  

  double RMAXVec[8];
  RMAXVec[0] = 70;
  RMAXVec[1] = 70;
  RMAXVec[2] = 70;
  RMAXVec[3] = 40;
  RMAXVec[4] = 40;
  RMAXVec[5] = 80;
  RMAXVec[6] = 80;
  RMAXVec[7] = 60; 

  double Z_Values[8];
  Z_Values[0] =-30;
  Z_Values[1] =-20;
  Z_Values[2] =-10;
  Z_Values[3] =  0;
  Z_Values[4] = 10;
  Z_Values[5] = 20;
  Z_Values[6] = 30;
  Z_Values[7] = 40;

  double Phi_Values[2];
   Phi_Values[0]=0*deg;

  Phi_Values[1]=45.*deg;
  // Phi_Values[0]=0;
  // Phi_Values[1]=2*pi;
  
  G4cout << "\n=======     Polyhedra test      ========";

  G4Polyhedra *MyPCone = new G4Polyhedra ("MyPCone",
                                                    Phi_Values[0],
						    Phi_Values[1],
						    5        ,
						    8        ,
						    Z_Values ,
						    RMINVec  ,
						    RMAXVec   );
  
  G4cout << "\n\nPolyhedra is  created ! "<<G4endl;
  // -> Check methods :
  //  - Inside
  //  - DistanceToIn
  //  - DistanceToOut

  EInside in;
  
  G4cout<<"\n\n==================================================";
  G4ThreeVector pt(0, -100, 24);
  G4int y;
  for (y = -100; y<=100; y+=10)
  {
    pt.setY(y);
    in = MyPCone->Inside(pt);
    
    G4cout << "\nx=" << pt.x() << "  y=" << pt.y() << "  z=" << pt.z();
    
    if( in == kInside )
      G4cout <<" is inside";
    else
      if( in == kOutside )
	G4cout <<" is outside";
      else
	G4cout <<" is on the surface";
  }

  G4cout<<"\n\n==================================================";
  G4ThreeVector start( 0, 0, -30);
  G4ThreeVector dir(1./std::sqrt(2.), 1./std::sqrt(2.), 0);
  G4double   d;
  G4int z;
  
  G4cout<<"\nPdep is (0, 0, z)";
  G4cout<<"\nDir is (1, 1, 0)\n";

  for(z=-30; z<=50; z+=5)
  {
    start.setZ(z);

    in = MyPCone->Inside(start);
    G4cout<< "x=" << start.x() << "  y=" << start.y() << "  z=" << start.z();
    
    if( in == kInside )
    {
      G4cout <<" is inside";

      d = MyPCone->DistanceToOut(start, dir);
      G4cout<<"  distance to out="<<d;
      d = MyPCone->DistanceToOut(start);
      G4cout<<"  closest distance to out="<<d<<G4endl;
    }
    else if( in == kOutside ) 
    {
      G4cout <<" is outside";

      d = MyPCone->DistanceToIn(start, dir);
      G4cout<<"  distance to in="<<d;
      d = MyPCone->DistanceToIn(start);
      G4cout<<"  closest distance to in="<<d<<G4endl;
    }
    else
      G4cout <<" is on the surface"<<G4endl;

  }

  G4cout<<"\n\n==================================================";
  G4ThreeVector start2( 0, -100, -30);
  G4ThreeVector dir2(0, 1, 0);
  G4double   d2;

  G4cout<<"\nPdep is (0, -100, z)";
  G4cout<<"\nDir is (0, 1, 0)\n";

  for(z=-30; z<=50; z+=5)
  {
    G4cout<<"  z="<<z;
    start2.setZ(z);
    d2 = MyPCone->DistanceToIn(start2, dir2);
    G4cout<<"  distance to in="<<d2;
    d2 = MyPCone->DistanceToIn(start2);
    G4cout<<"  distance to in="<<d2<<G4endl;
  }

  G4cout<<"\n\n==================================================";
  G4ThreeVector start3( 0, 0, -50);
  G4ThreeVector dir3(0, 0, 1);
  G4double   d3;

  G4cout<<"\nPdep is (0, y, -50)";
  G4cout<<"\nDir is (0, 0, 1)\n";

  for(y=-0; y<=90; y+=5)
  {
    G4cout<<"  y="<<y;
    start3.setY(y);
    d3 = MyPCone->DistanceToIn(start3, dir3);
    G4cout<<"  distance to in="<<d3<<G4endl;
  }
  //
  // Add checks in Phi direction
  // Point move in Phi direction for differents Z
  //
 G4cout<<"\n\n==================================================";
   
 for(z=-10; z<=50; z+=5)
{
  G4cout<<"\n\n===================Z="<<z<<"==============================";
     
  G4ThreeVector start4( 0, 0, z);
  G4double phi=45.*pi/180.*rad;
  G4ThreeVector dir4(std::cos(phi), std::sin(phi), 0);
  G4double   d4;

  G4cout<<"\nPdep is (0<<R<<50, phi, z)";
  G4cout<<"\nDir is (std::cos(phi), std::sin(phi), 0)\n";
  G4cout<<"Ndirection is="<<dir4 <<G4endl;

  for(y=-0; y<=50; y+=5)
  {
    
    start4.setX(y*std::cos(phi));
    start4.setY(y*std::sin(phi));
    G4cout<<"  R="<<y<<" with Start"<<start4;
    in = MyPCone->Inside(start4);
    if( in == kInside )
      {
       G4cout <<" is inside";
        d4 = MyPCone->DistanceToOut(start4, dir4);
         G4cout<<"  distance to out="<<d4;
         d4 = MyPCone->DistanceToOut(start4);
         G4cout<<" closest distance to out="<<d4<<G4endl;
	}
    else
      if( in == kOutside )
	{
         G4cout <<" is outside";
          d4 = MyPCone->DistanceToIn(start4, dir4);
         G4cout<<"  distance to in="<<d4;
         d4 = MyPCone->DistanceToIn(start4);
         G4cout<<" closest distance to in="<<d4<<G4endl;
	}
      else
	{G4cout <<" is on the surface";
         d4 = MyPCone->DistanceToIn(start4, dir4);
         G4cout<<"  distance to in="<<d4;
         d4 = MyPCone->DistanceToIn(start4);
         G4cout<<" closest distance to in="<<d4<<G4endl;
	}
    
  }
}
 //
 // Add checks in Phi direction
 // Point move in X direction for differents Z
 // and 'schoot' on rhi edge
 G4cout<<"\n\n==================================================";
 
 for(z=-10; z<=50; z+=5)
 {
  G4cout<<"\n\n===================Z="<<z<<"==============================";
     
  G4ThreeVector start5( 0., 1, z);
  G4ThreeVector dir5(0,-1, 0);
  G4double   d5;

  G4cout<<"\nPdep is (0<<X<<50, 1, z)";
  G4cout<<"\nDir is (0, -1, 0)\n";
  G4cout<<"Ndirection is="<<dir5 <<G4endl;

  for(y=-0; y<=50; y+=5)
  {
    
    start5.setX(y);
    G4cout<<" Start"<<start5;
    in = MyPCone->Inside(start5);
    if( in == kInside )
      {
       G4cout <<" is inside";
       d5 = MyPCone->DistanceToOut(start5, dir5);
       G4cout<<"  distance to out="<<d5;
       d5 = MyPCone->DistanceToOut(start5);
       G4cout<<" closest distance to out="<<d5<<G4endl;
      }
    else
      if( in == kOutside )
        {
	 G4cout <<" is outside";
         d5 = MyPCone->DistanceToIn(start5, dir5);
         G4cout<<"  distance to in="<<d5;
         d5 = MyPCone->DistanceToIn(start5);
         G4cout<<" closest distance to in="<<d5<<G4endl;
        }
      else
        {
	 G4cout <<" is on the surface";
         d5 = MyPCone->DistanceToIn(start5, dir5);
         G4cout<<"  distance to in="<<d5;
         d5 = MyPCone->DistanceToIn(start5);
         G4cout<<" closest distance to in="<<d5<<G4endl;
        }
    
  }
 }


 // Asserts
 G4ThreeVector p1,p2,p3,p4,p5,p6,dirx,diry,dirz;
 p1=G4ThreeVector(0,0,-5); 
 p2=G4ThreeVector(50,0,40);
 p3=G4ThreeVector(5,1,20 ); 
 p4=G4ThreeVector(45,5,30);
 p5=G4ThreeVector(0,0,30); 
 p6=G4ThreeVector(41,0,10);

 dirx=G4ThreeVector(1,0,0);
 diry=G4ThreeVector(0,1,0);
 dirz=G4ThreeVector(0,0,1);

 // Inside
  assert(MyPCone->Inside(p1) == kSurface);
  assert(MyPCone->Inside(p2) == kSurface);
  assert(MyPCone->Inside(p3) == kInside);
  assert(MyPCone->Inside(p4) == kInside);
  assert(MyPCone->Inside(p5) == kOutside);
  assert(MyPCone->Inside(p6) == kOutside);
 // DistanceToIn
 
  assert(std::fabs((MyPCone->DistanceToIn(p1,dirx))) < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToIn(p1,-diry)))< kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToIn(p2,diry))) < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToIn(p5,dirx)  -40.12368793931)) < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToIn(p6,-dirx)  -0.87631206069)) < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToIn(p6,dirz)   -0.218402670765))< kCarTolerance);
 // DistanceToOut
  assert(std::fabs((MyPCone->DistanceToOut (p1,-dirx)))               < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToOut (p3,-diry) -1.) )          < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToOut (p3,dirz)  -1.27382374146))< kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToOut(p4,dirz)  -10.))           < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToOut(p4,dirx)  -34.8538673445)) < kCarTolerance);
  assert(std::fabs((MyPCone->DistanceToOut(p4,diry)  -40.))           < kCarTolerance);



  return EXIT_SUCCESS;
}

