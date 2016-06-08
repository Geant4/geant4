#ifndef G4UserLimitsForRD_h
#define G4UserLimitsForRD_h 1

#include "globals.hh"
#include "G4UserLimits.hh"

class G4Track;

////////////////////////////////////////////////////////////////////////////////
//
class G4UserLimitsForRD : public G4UserLimits
{
  // class description 
  //
  //   Derived class from G4UserLimits.  It is specific for
  //   RDM to select and deselect volumes in which RDM is applied
  //
  // class  description - end

public:
    G4UserLimitsForRD ( const G4String aType="RD") :
      G4UserLimits (aType)
  {  
#ifdef G4VERBOSE
    if (GetVerboseLevel()>1)
      G4cerr <<"G4UserLimitsForRD constructor" <<endl;
#endif
  }
  ~G4UserLimitsForRD () {;}

  //  G4bool IsRDActive(const G4Track&);
};
#endif

