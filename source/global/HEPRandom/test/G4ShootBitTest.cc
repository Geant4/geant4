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
// $Id: G4ShootBitTest.cc,v 1.7 2006-06-29 19:00:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

// Tests RandFlat::fireBit() member function
// Peter Urban, 5th Sep 1996

#include "G4ios.hh"
#include <iomanip>
#include "Randomize.hh"
#include "G4Timer.hh"

int main()
{
  HepRandomEngine* e= HepRandom::getTheEngine();
  RandFlat r(*e);
  
  const long NumTest= 10000;  

  G4cout << "Printed value: (relative frequency - probability)" << G4endl;
  G4cout << NumTest << " random numbers" << G4endl << G4endl;

  long sum= 0;
  long i, j;

  for (i=0; i<NumTest; i++) {
    sum+= r.fireBit();
  }
  G4cout << "fireBit                   " << double(sum)/NumTest-0.5 << G4endl;  
  sum= 0;
  for (i=0; i<NumTest; i++) {
    sum+= r.shootBit();
  }
  G4cout << "shootBit                  " << double(sum)/NumTest-0.5 << G4endl;  
  sum= 0;
  for (i=0; i<NumTest; i++) {
    sum+= r.shootBit(e);
  }
  G4cout << "shootBit(engine)          " << double(sum)/NumTest-0.5 << G4endl;  

 
  const long Num= 1000000;
  int iii;  long lll;

  G4Timer theTimer;
  G4cout << G4endl << "Performance test: " << G4endl;

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
  G4cout << theTimer.GetUserElapsed() << G4endl << G4endl;
  
  return 0;
}



