#include "G4VStatMFMacroCluster.hh"



// Copy constructor
G4VStatMFMacroCluster::G4VStatMFMacroCluster(const G4VStatMFMacroCluster & right)
{
	G4Exception("G4VStatMFMacroCluster::copy_constructor meant to not be accessable");
}

// Operators

G4VStatMFMacroCluster & G4VStatMFMacroCluster::
operator=(const G4VStatMFMacroCluster & right)
{
	G4Exception("G4VStatMFMacroCluster::operator= meant to not be accessable");
	return *this;
}


G4bool G4VStatMFMacroCluster::operator==(const G4VStatMFMacroCluster & right) const
{
//	G4Exception("G4VStatMFMacroCluster::operator== meant to not be accessable");
	return false;
}
 

G4bool G4VStatMFMacroCluster::operator!=(const G4VStatMFMacroCluster & right) const
{
//	G4Exception("G4VStatMFMacroCluster::operator!= meant to not be accessable");
	return true;
}


G4double G4VStatMFMacroCluster::CalcInvLevelDensity(void)
{
	// Calculate Inverse Density Level
	// Epsilon0*(1 + 3 /(Af - 1))
	if (theA == 1) return 0.0;
	else return
		G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(G4double(theA) - 1.0));
}

