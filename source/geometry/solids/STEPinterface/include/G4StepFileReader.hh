// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StepFileReader.hh,v 1.2 1999-12-15 14:50:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef G4STEPFILEREADER_HH
#define G4STEPFILEREADER_HH

#include "globals.hh"
#include "instmgr.h"

class G4StepFileReader
{
public:
  virtual void ReadSTEPFile(G4String)=0;
  virtual void SaveSTEPFile()=0;
  virtual void UpdateSTEPFile()=0;
  virtual InstMgr GetInstanceManager()=0;
};

#endif
