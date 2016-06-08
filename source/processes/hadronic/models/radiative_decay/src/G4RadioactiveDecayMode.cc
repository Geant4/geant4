
#include "G4RadioactiveDecayMode.hh"

G4std::istream &operator >> (G4std::istream &s, G4RadioactiveDecayMode &q)
{
  G4String a;
  s >> a;
  if (a == "IT")
    {q = IT;}
  else if (a == "BetaMinus")
    {q = BetaMinus;}
  else if (a == "BetaPlus")
    {q = BetaPlus;}
  else if (a == "KshellEC")
    {q = KshellEC;}
  else if (a == "LshellEC")
    {q = LshellEC;}
  else if (a == "MshellEC")
    {q = MshellEC;}
  else if (a == "Alpha")
    {q = Alpha;}
  else
    {q = ERROR;}
  return s;
}
////////////////////////////////////////////////////////////////////////////////

