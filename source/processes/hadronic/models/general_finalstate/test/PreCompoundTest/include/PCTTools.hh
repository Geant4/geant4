#ifndef PCTTools_hh
#define PCTTools_hh 1

#include "globals.hh"
#include <strstream>

inline G4double StrToDouble(const G4String& s)
{
  std::istringstream tonum(s.c_str());
  G4double result;
  tonum >> result;
  return result;
}

inline G4int StrToInt(const G4String& s)
{
  std::istringstream tonum(s.c_str());
  G4int result;
  tonum >> result;
  return result;
}

inline G4String DoubleToStr(const G4double n)
{
    std::ostringstream conv;
    conv << n;
    return conv.str();
}


inline G4String IntToStr(const G4int n)
{
    std::ostringstream conv;
    conv << n;
    return conv.str();
}

#endif // PCTTools_hh
