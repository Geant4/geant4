//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: OrderedTableTest.cc,v 1.6 2001-09-17 08:17:59 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ----------------------------------------------------------------------
//
// This program shows how to use the G4OrderedTable. 
//

#include "G4ios.hh"
#include "G4OrderedTable.hh"

int main()
{
  size_t I, J;
  const size_t Imax=10;
  G4cout.precision(3); 

 //
 // Create a G4OrderedTable object
 // 

  G4OrderedTable aTable(Imax);
  G4OrderedTable::iterator pl = aTable.begin();

  for(I=0; I<Imax; I++)
  {
    G4DataVector* aVector = new G4DataVector();
    *pl++ = aVector;
    for(J=0; J<=I; J++)
      aVector->push_back(G4double(J));
  }

 // Now access the data contained in the table

  for (I=0; I<Imax; I++)
  {
    G4cout << G4endl << G4endl << " I= " << I << "  Data= ";
    for(J=0; J<=I; J++)
      G4cout << (*aTable[I])[J] << " ";
  }


 // Store in file in ascii mode 
  aTable.Store("OrdTable.asc",true) ; 

 // clear Ordered Table
  aTable.clear();

 // Retrieve from file
  aTable.Retrieve("OrdTable.asc",true) ; 

 // Print Out 
  G4cout <<   G4endl << G4endl;
  G4cout << aTable ;

 //
 // Create a G4OrderedTable object by pointer 
 //

  G4OrderedTable*  aTablePtr = new G4OrderedTable(Imax); 
  pl = aTablePtr->begin();

  for(I=0; I<Imax; I++)
  {
    G4DataVector* aVector = new G4DataVector(I+1);
    *pl++ = aVector;
    for(J=0; J<=I; J++)
      (*aVector)[J]= G4double(J);
  }

 // Now access the data contained in the table 

  for (I=0; I<Imax; I++)
  {
    G4cout << G4endl << G4endl << " I= " << I << "   Data= ";
    for(J=0; J<=I; J++)
      G4cout << (*(*aTablePtr)[I])[J] << " ";
  }
  G4cout << G4endl;

  // Store in file in binary mode 
  aTablePtr->Store("OrdTable.dat") ; 

  aTablePtr->clearAndDestroy();

 // Retrieve from file
  aTablePtr->Retrieve("OrdTable.dat") ; 

 // Print Out 
  G4cout <<   G4endl << G4endl;
  G4cout << *aTablePtr ;


  aTable.clearAndDestroy();
  aTablePtr->clearAndDestroy();

  return EXIT_SUCCESS;
}
