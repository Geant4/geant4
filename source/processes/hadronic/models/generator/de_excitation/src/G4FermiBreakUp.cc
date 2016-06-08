// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4FermiBreakUp.hh"


G4FermiBreakUp::G4FermiBreakUp()
{
}

G4FermiBreakUp::G4FermiBreakUp(const G4FermiBreakUp &right)
{
  G4Exception("G4FermiBreakUp::copy_constructor meant to not be accessable");
}


G4FermiBreakUp::~G4FermiBreakUp()
{
}


const G4FermiBreakUp & G4FermiBreakUp::operator=(const G4FermiBreakUp &right)
{
  G4Exception("G4FermiBreakUp::operator= meant to not be accessable");
  return *this;
}


G4bool G4FermiBreakUp::operator==(const G4FermiBreakUp &right) const
{
  return false;
}

G4bool G4FermiBreakUp::operator!=(const G4FermiBreakUp &right) const
{
  return true;
}



G4FragmentVector * G4FermiBreakUp::BreakItUp(const G4Fragment &theNucleus)
{
	// CHECK that Excitation Energy > 0
	if (theNucleus.GetExcitationEnergy() <= theNucleus.GetBindingEnergy()) {
		G4FragmentVector * theResult = new G4FragmentVector;
		theResult->insert(new G4Fragment(theNucleus));
		return theResult;
	}

	// Total energy of nucleus in nucleus rest frame 
	G4double TotalEnergyRF = theNucleus.GetExcitationEnergy() +
		G4ParticleTable::GetParticleTable()->GetIonTable()->
		GetIonMass(theNucleus.GetZ(),theNucleus.GetA());

	G4FermiConfigurationList theConfigurationList;


	// Split the nucleus
	G4bool Split = theConfigurationList.Initialize(theNucleus.GetA(), 
																  theNucleus.GetZ(),
																  TotalEnergyRF);
	if ( !Split ) {
		G4FragmentVector * theResult = new G4FragmentVector;
		theResult->insert(new G4Fragment(theNucleus));

		return theResult;
	}

	// Chose a configuration
	G4FermiConfiguration theConfiguration(theConfigurationList.ChooseConfiguration());
  

	// Get the fragments corresponding to chosen configuration.
	G4FragmentVector * theResult = theConfiguration.GetFragments(theNucleus);

	return theResult;
	
}


