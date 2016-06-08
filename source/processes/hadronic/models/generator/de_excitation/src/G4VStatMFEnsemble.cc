#include "G4VStatMFEnsemble.hh"



// Copy constructor
G4VStatMFEnsemble::G4VStatMFEnsemble(const G4VStatMFEnsemble & right)
{
	G4Exception("G4VStatMFEnsemble::copy_constructor meant to not be accessable");
}

// Operators

G4VStatMFEnsemble & G4VStatMFEnsemble::
operator=(const G4VStatMFEnsemble & right)
{
	G4Exception("G4VStatMFEnsemble::operator= meant to not be accessable");
	return *this;
}


G4bool G4VStatMFEnsemble::operator==(const G4VStatMFEnsemble & right) const
{
	G4Exception("G4VStatMFEnsemble::operator== meant to not be accessable");
	return false;
}
 

G4bool G4VStatMFEnsemble::operator!=(const G4VStatMFEnsemble & right) const
{
	G4Exception("G4VStatMFEnsemble::operator!= meant to not be accessable");
	return true;
}

