// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ShootBitTest.cc,v 1.2 1999-11-16 17:31:42 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Tests RandFlat::fireBit() member function
// Peter Urban, 5th Sep 1996

#include "G4ios.hh"
#include <iomanip.h>
#include "Randomize.hh"
#include "G4Timer.hh"

int main()
{
  HepRandomEngine* e= HepRandom::getTheEngine();
  RandFlat r(*e);
  
  const long NumTest= 10000;  

  G4cout << "Printed value: (relative frequency - probability)" << endl;
  G4cout << NumTest << " random numbers" << endl << endl;

  long sum= 0;
  long i, j;

  for (i=0; i<NumTest; i++) {
    sum+= r.fireBit();
  }
  G4cout << "fireBit                   " << double(sum)/NumTest-0.5 << endl;  
  sum= 0;
  for (i=0; i<NumTest; i++) {
    sum+= r.shootBit();
  }
  G4cout << "shootBit                  " << double(sum)/NumTest-0.5 << endl;  
  sum= 0;
  for (i=0; i<NumTest; i++) {
    sum+= r.shootBit(e);
  }
  G4cout << "shootBit(engine)          " << double(sum)/NumTest-0.5 << endl;  

 
  const long Num= 1000000;
  int iii;  long lll;

  G4Timer theTimer;
  G4cout << endl << "Performance test: " << endl;

  G4cout << "shootBit() vs G4UniformRand()<0.5: ";
  theTimer.Start();
  for (i=0; i<Num; i++)
    iii = RandFlat::shootBit();
  theTimer.Stop();
  G4cout << theTimer.GetUserElapsed() << " vs ";
  theTimer.Start();
  for (i=0; i<Num; i++)
    G4UniformRand()<0.5;
  theTimer.Stop();
  G4cout << theTimer.GetUserElapsed() << endl << endl;
  
  return 0;
}



