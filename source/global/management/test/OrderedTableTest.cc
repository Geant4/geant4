// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: OrderedTableTest.cc,v 1.3 1999-11-23 15:00:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//
// This program shows how to use the G4OrderedTable. 
//

#include "G4ios.hh"
#include "G4OrderedTable.hh"

G4int main() {

  G4int Imax=10;
  G4cout.precision(3); 

 //
 // Create a G4OrderedTable object
 // 

  G4OrderedTable  aTable(Imax); 

  for(G4int I=0; I<Imax; I++){
    aTable(I) = new G4ValVector(I+1);    // aTable.insertAt(I,new G4ValVector(I+1));
    for(G4int J=0; J<=I; J++) {
      (*aTable(I))(J) = G4double(J);     // aTable(I)->insertAt(J,G4double(J));
    }
  }

 // Now access the data contained in the table

  for (I=0; I<Imax; I++){
    G4cout << G4endl << G4endl << " I= " << I << "  Data= ";
    for(G4int J=0; J<=I; J++) {
      G4cout << (*aTable(I))(J) << " ";
    }
  }

 //
 // Create a G4OrderedTable object by pointer 
 //

  G4OrderedTable*  aTablePtr = new G4OrderedTable(Imax); 

  for(I=0; I<Imax; I++){
    (*aTablePtr)(I) = new G4ValVector(I+1);  // aTablePtr->insertAt(I,new G4ValVector(I+1))
    for(G4int J=0; J<=I; J++) {
      (*(*aTablePtr)(I))(J) = G4double(J);   // (*aTablePtr)(I)->insertAt(J,G4double(J))
    }
  }

 // Now access the data contained in the table 

  for (I=0; I<Imax; I++){
    G4cout << G4endl << G4endl << " I= " << I << "   Data= ";
    for(G4int J=0; J<=I; J++) {
      G4cout << (*(*aTablePtr)(I))(J) << " ";
    }
  }
  G4cout << G4endl;

  return EXIT_SUCCESS;
}

