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

// Other Class Headers
#include "globals.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"

// C++ Headers

//RW Headers
#include <rw/tpslist.h>

// Class Declarations
class G4VDataFile;

class G4EpdlTables{

public:
  
  // Constructors
  G4EpdlTables(G4VDataFile& DFile);
  
  // Destructor 
  ~G4EpdlTables();

  // Member Functions

  // search the data table in the file
  void FillDataTable();
  void FillTheTable(G4int nemEl = 0);

  // inline member functions

  inline G4PhysicsTable* GetFstDataTable(){return new G4PhysicsTable(*theDataTable1);};
  inline G4PhysicsTable* GetSndDataTable(){return theDataTable2;};
  inline G4PhysicsTable* GetTrdDataTable(){return theDataTable3;};

  RWTPtrSlist< RWTPtrSlist<G4DataVector> >* GetGlobalList();

  void MakeNtuple(const char* fname);

private:

protected:


private:

  G4PhysicsTable* theDataTable1;
  G4PhysicsTable* theDataTable2;
  G4PhysicsTable* theDataTable3;
  RWTPtrSlist< RWTPtrSlist<G4DataVector> > allElementList;

  G4VDataFile& datfile;

};

#endif // G4EpdlTables









