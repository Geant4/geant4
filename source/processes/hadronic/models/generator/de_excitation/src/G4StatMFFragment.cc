
#include "G4StatMFFragment.hh"



// Copy constructor

G4StatMFFragment::G4StatMFFragment(const G4StatMFFragment & right)
{
  G4Exception("G4StatMFFragment::copy_constructor menat to not be accessable");
}


// Operators

const G4StatMFFragment & G4StatMFFragment::operator=(const G4StatMFFragment & right)
{
  G4Exception("G4StatMFFragment::operator= meant to not be accessable");
  return *this;
}


G4bool G4StatMFFragment::operator==(const G4StatMFFragment & right) const 
{
  return false;
}
 

G4bool G4StatMFFragment::operator!=(const G4StatMFFragment & right) const 
{
  return true;
}

