#ifndef G4Fuzzy_hh
#define G4Fuzzy_hh

#include <utility>

class G4Fuzzy : public std::pair<bool, bool>
{
  public:
  G4Fuzzy() 
  {
    first=true;
    second=false;
  }
  G4Fuzzy(bool aB) 
  {
    first=false;
    second=aB; 
  }
};

#endif
