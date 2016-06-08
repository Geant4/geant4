// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4EpdlTables
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 2 February 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------
#ifndef G4EpdlTables_hh
#define G4EpdlTables_hh

// Base Class Header
#include "G4VTables.hh"

// Other Class Headers
#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4Data.hh"
#include "G4FirstLevel.hh"
#include "G4SecondLevel.hh"

// C++ Headers

// Class Declarations
class G4VDataFile;

class G4EpdlTables: public G4VTables{

public:
  
  // Constructors
  G4EpdlTables(G4VDataFile& DFile);
  
  // Destructor 
  ~G4EpdlTables();

  // Member Functions

  // search the data table in the file
  void FillDataTable();
  G4SecondLevel* FillTheTable(G4int nemEl = 0);

  // inline member functions

  inline G4PhysicsTable* GetFstDataTable(){return theDataTable1;};
  inline G4PhysicsTable* GetSndDataTable(){return theDataTable2;};
  inline G4PhysicsTable* GetTrdDataTable(){return theDataTable3;};

  //G4SecondLevel* GetGlobalList();

protected:


private:

  G4PhysicsTable* theDataTable1;
  G4PhysicsTable* theDataTable2;
  G4PhysicsTable* theDataTable3;
  //  G4SecondLevel* allElementList;

  G4VDataFile& datfile;

};

#endif // G4EpdlTables









