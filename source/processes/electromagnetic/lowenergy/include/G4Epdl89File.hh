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
//      File name:     G4Epdl89File
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 2 February 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------
#ifndef G4Epdl89File_hh
#define G4Epdl89File_hh

// Other Class Headers
#include "globals.hh"
#include "G4VDataFile.hh"

// C++ Headers
#include <fstream.h>

// Class Declarations

class G4Epdl89File: public G4VDataFile{

public:
  
  // Constructors
  G4Epdl89File(const G4String&, G4int*);
  
  // Destructor 
  ~G4Epdl89File();

  virtual G4bool FindTheProcess(); 
  virtual G4bool FindTheElement(G4int);
  virtual G4bool FindOneElemProc(G4int&);

protected:

  // Member Functions

  void GetDataValues(G4DataVector& valList);

  G4double GetOneData(const char*);

private:

  G4int* _flags;                     // data flags
};

#endif // G4Epdl89File





