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
#include "G4BREPSolidSphere.hh"
#include "G4Timer.hh"



int main(int argc, char **argv)
{
  G4Timer timer;

  G4ThreeVector tStart(19000,0,10.1);
  G4ThreeVector tDir(-1,0,0);
  G4ThreeVector tStart2(0,0,10);
  G4ThreeVector pt(0,0,50);
  G4ThreeVector pt2(100000,0,50);
  double Dist;


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
  G4cout << "\nStart 0,0,10";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MySphere->DistanceToOut(pt);  
  G4cout << "\nPoint 0,0,50";
  G4cout << "\nDist to out: " << Dist ;

  Dist = MySphere->DistanceToIn(pt2);  
  G4cout << "\nPoint 100000,0,50";
  G4cout << "\nDist to in: " << Dist << "\n";

  G4cout << endl << endl;

  return EXIT_SUCCESS;
}

