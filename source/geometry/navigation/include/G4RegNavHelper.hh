#ifndef G4RegNavHelper_HH
#define G4RegNavHelper_HH

#include <vector>
#include "globals.hh"

class G4RegNavHelper
{
public:
  G4RegNavHelper(){};
  ~G4RegNavHelper(){};
  
  static void ClearStepLengths();
  static void AddStepLength(G4int copyNo, G4double slen );

  static std::vector< std::pair<G4int,G4double> > theStepLengths;

};

#endif
