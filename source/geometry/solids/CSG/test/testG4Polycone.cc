// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// First Polycone test,             created by J. Apostolakis,  12 Feb 99
// modification of BREP PCon test,  created by L. Broglia,      20 Oct 98
//      which was derived from G4Gerep test by J. Sulkimo
//
// Cvs version: $ Id $
// Cvs tag :  $ Name $

#include <stdio.h>
#include <math.h>
#include <fstream.h>
#include <stdlib.h>
#include "G4ios.hh" 
// #include "G4BREPSolid.hh"
#include "G4Polycone.hh"
#include "G4Timer.hh"

#include "G4Vector3D.hh"
#include "G4Point3D.hh"

int main(int argc, char **argv)
{
  G4Timer timer;
  
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
  Z_Values[0] =-20;
  Z_Values[1] =-10;
  Z_Values[2] =-10;
  Z_Values[3] =  0;
  Z_Values[4] = 10;
  Z_Values[5] = 20;
  Z_Values[6] = 30;
  Z_Values[7] = 40;
  
  G4cout << "\n=======     Polycone test      ========";

  G4Polycone *MyPCone = new G4Polycone ("MyPCone",
						    0        ,
						    2*pi     ,
						    8        ,
						    Z_Values ,
						    RMINVec  ,
						    RMAXVec   );
  
  G4cout << "\n\nPCone created ! "<<endl;
  // -> Check methods :
  //  - Inside
  //  - DistanceToIn
  //  - DistanceToOut

  
  EInside in;
  
  G4cout<<"\n\n==================================================";
  G4Point3D  pt(0, -100, 24);
  for (G4int y = -100; y<=100; y+=10)
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
  G4Point3D  start( 0, 0, -30);
  G4Vector3D dir(1, 1, 0);
  G4double   d;
  
  G4cout<<"\nPdep is (0, 0, z)";
  G4cout<<"\nDir is (1, 1, 0)\n";

  for(G4double z=-30; z<=50; z+=5)
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
      G4cout<<"  closest distance to out="<<d<<endl;
    }
    else if( in == kOutside ) 
    {
      G4cout <<" is outside";

      d = MyPCone->DistanceToIn(start, dir);
      G4cout<<"  distance to in="<<d;
      d = MyPCone->DistanceToIn(start);
      G4cout<<"  closest distance to in="<<d<<endl;
    }
    else
      G4cout <<" is on the surface"<<endl;

  }

  G4cout<<"\n\n==================================================";
  G4Point3D  start2( 0, -100, -30);
  G4Vector3D dir2(0, 1, 0);
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
    G4cout<<"  distance to in="<<d2<<endl;
  }

  G4cout<<"\n\n==================================================";
  G4Point3D  start3( 0, 0, -50);
  G4Vector3D dir3(0, 0, 1);
  G4double   d3;

  G4cout<<"\nPdep is (0, y, -50)";
  G4cout<<"\nDir is (0, 0, 1)\n";

  for(y=-0; y<=90; y+=5)
  {
    G4cout<<"  y="<<y;
    start3.setY(y);
    d3 = MyPCone->DistanceToIn(start3, dir3);
    G4cout<<"  distance to in="<<d3<<endl;
  }
  
  
  return EXIT_SUCCESS;
}

