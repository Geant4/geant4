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


#include "globals.hh"
#include "G4BREPSolid.hh"
#include "G4BREPSolidBox.hh"
#include "G4Timer.hh"


int main()
{
  G4Timer timer;

  G4ThreeVector tStart(19000,0,10.1);
  G4ThreeVector tDir(-1,0,0);
  G4ThreeVector tStart2(0,0,10);
  G4ThreeVector pt(0,0,50);
  G4ThreeVector pt2(100000,0,50);
  double d1, d2, d3, d4 ;

cout << "\n ============   Box test   ================"; 

  G4BREPSolidBox *myCalBox = new G4BREPSolidBox 
    ( "MyBox",
      G4Point3D(-1500, -1500, -1000),
      G4Point3D(-1500, -1500,  1000),       
      G4Point3D(-1500,  1500,  1000),
      G4Point3D(-1500,  1500, -1000),
      G4Point3D( 1500, -1500, -1000),       
      G4Point3D( 1500, -1500,  1000),
      G4Point3D( 1500,  1500,  1000),
      G4Point3D( 1500,  1500, -1000) );

  cout << "\n\nBox created ! ";
  
  // -> Check methods :
  //  - Inside
  //  - DistanceToIn
  //  - DistanceToOut

  G4Point3D Pt[4];
  Pt[0] = G4Point3D(    1,    1,  100);
  Pt[1] = G4Point3D( 1000, 1000, 5000);
  Pt[2] = G4Point3D(-1000,-1000,-5000);
  Pt[3] = G4Point3D(    0,    0, -100);

  G4Vector3D Dir[4];
  Dir[0] = G4Vector3D(    1,    0,    0);
  Dir[1] = G4Vector3D(    0,    1,    0);
  Dir[2] = G4Vector3D(    0,    0,    1);
  Dir[3] = G4Vector3D(    1,    1,   -1);

  EInside in[4];
  G4double dist[4][2];

  G4cout<<"\n\n";

  for (G4int a = 0; a<4; a++)
  {
    in[a] = myCalBox->Inside(Pt[a]);

    G4cout<<"----------------------------\n\n";

    G4cout<<"x="<<Pt[a].x()
	  <<"  y="<<Pt[a].y()<<"  z="<<Pt[a].z();
      
    if( in[a] == kInside )
    {
      G4cout <<" is inside"<<endl;

      dist[a][1] = myCalBox->DistanceToOut(Pt[a]);
      cout<<"\nDistance to out is :"<<dist[a][1]<<endl;

      cout << "\nDir   : x=" << Dir[a].x() 
	   << " y=" << Dir[a].y() 
	   << " z=" << Dir[a].z()<<endl;
      dist[a][2] = myCalBox->DistanceToOut(Pt[a], Dir[a]);
      cout<<"Distance to out is :"<<dist[a][2]<<endl;

    }
    else
    {
      G4cout <<" is outside"<<endl;

      dist[a][1] = myCalBox->DistanceToIn(Pt[a]);
      cout<<"\nDistance to in is :"<<dist[a][1]<<endl;

      cout << "\nDir   : x=" << Dir[a].x() 
	   << " y=" << Dir[a].y() 
	   << " z=" << Dir[a].z()<<endl;
      dist[a][2] = myCalBox->DistanceToIn(Pt[a], Dir[a]);
      cout<<"Distance to in is :"<<dist[a][2];
    }
    cout<<endl;
  }

  return EXIT_SUCCESS;
}

