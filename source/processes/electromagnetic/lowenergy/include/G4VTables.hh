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
//      File name:     G4VTables
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 2 February 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------
#ifndef G4VTables_hh
#define G4VTables_hh

#include "globals.hh"

class G4VTables{

public:
  
  // Constructors
  G4VTables();
  
  // Destructor 
  ~G4VTables();

  // Member Functions

  // search the data table in the file
  virtual void FillDataTable() = 0;
  virtual void FillTheTable(G4int nemEl = 0) = 0;

protected:


private:

};

#endif // G4VTables









