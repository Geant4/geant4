#ifndef __G4STEPENTITY
#define __G4STEPENTITY
#include "globals.hh"
#include "G4OrderedTable.hh"

class G4STEPEntity
{
public:
  virtual G4String GetEntityType()=0;
};



#endif
