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
//      File name:     G4VDataFile
//
//      Author:        Alessandra Forti (Alessandra.Forti@cern.ch)
// 
//      Creation date: 2 February 1999
//
//      Modifications: 
//      
// -------------------------------------------------------------------
#ifndef G4VDataFile_hh
#define G4VDataFile_hh

// Other Class Headers
#include "globals.hh"
#include "G4DataVector.hh"

// C++ Headers
#include <fstream.h>

class G4VDataFile{

public:
  
  // Constructors
  G4VDataFile(const G4String& file);
  
  // Destructor 
  ~G4VDataFile();

  // Member Functions
  void OpenFile();
  void CloseFile();
  void Eof();
  G4bool IsOpen();
  void SeekPos(streampos);
  streampos TellPos();

  void GetLine();
  G4int LineLength();
  char* GetBuf();

  virtual G4bool FindTheProcess() = 0;
  virtual G4bool FindTheElement(G4int) = 0;
  virtual G4bool FindOneElemProc(G4int& subsh) = 0;
  virtual void GetDataValues(G4DataVector& valList) = 0;

protected:

  // Member Functions
  void SetBufferSize(G4int);

private:

  const G4String& _filename;             // data file name
  ifstream _istr;
  G4int _bufSize;
  char* buf;
};

#endif // G4VDataFile





















