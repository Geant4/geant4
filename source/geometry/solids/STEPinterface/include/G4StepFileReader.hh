// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4StepFileReader.hh,v 1.3 2000-01-21 13:45:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// Class G4StepFileReader
//
// Class description:
//
//

// Authors: J.Sulkimo, P.Urban.
// Revisions by: L.Broglia, G.Cosmo.
//
// History:
//   18-Nov-1999: First step of re-engineering - G.Cosmo
// ----------------------------------------------------------------------
#ifndef G4STEPFILEREADER_HH
#define G4STEPFILEREADER_HH

#include <schema.h>
#include "globals.hh"

class G4StepFileReader
{
  public:
    virtual void ReadSTEPFile(G4String)=0;
    virtual void SaveSTEPFile()=0;
    virtual void UpdateSTEPFile()=0;
    virtual InstMgr GetInstanceManager()=0;
};

#endif
