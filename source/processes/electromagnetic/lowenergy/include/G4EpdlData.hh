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
//      File name:     G4EpdlData
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 29 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------
#ifndef G4EPDLDATA_HH
#define G4EPDLDATA_HH

// Other classes headers
#include "CLHEP/String/Strings.h"
#include "G4PhysicsTable.hh"
#include "globals.hh"

// C++ headers

// class declarations
class HepTupleManager;

class G4EpdlData{

public:
  
  // Constructors
  G4EpdlData(HepString & , G4double*);
  
  // Destructor 
  ~G4EpdlData();

  // member functions

  // search the data table in the file
  void FillDataTable(G4double unit1, G4double unit2);

  // Make ntuples of the data
  void MakeNtuple(const char*);

  // inline member functions
  inline G4PhysicsTable* GetFstDataTable(){return theDataTable1;};
  inline G4PhysicsTable* GetSndDataTable(){return theDataTable2;};
  inline G4PhysicsTable* GetTrdDataTable(){return theDataTable3;};

protected:

private:

  G4PhysicsTable* theDataTable1;
  G4PhysicsTable* theDataTable2;
  G4PhysicsTable* theDataTable3;
  HepString _filename;                      // data file name
  G4double*_userflags;                     // data flags

  // member functions
  void CharToInt(char* linebuf, G4double* local); // convert a char linebuf to a string 
  void GetDataValues(char*, G4double*);
};

#endif // G4EpdlData


















