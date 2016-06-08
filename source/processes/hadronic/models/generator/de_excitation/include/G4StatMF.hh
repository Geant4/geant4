#ifndef G4StatMF_h
#define G4StatMF_h 1

#include "g4rw/tvordvec.h"

#include "globals.hh"
#include "G4VMultiFragmentation.hh"
#include "G4VStatMFEnsemble.hh"
#include "G4StatMFMicroCanonical.hh"
#include "G4StatMFMacroCanonical.hh"
#include "G4StatMFChannel.hh"
#include "G4Fragment.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"



class G4StatMF : public G4VMultiFragmentation
{
public:
	// Default constructor
	G4StatMF();
	// Destructor
	~G4StatMF();

private:
	// Copy constructor	
	G4StatMF(const G4StatMF & right);

	// Operators
	G4StatMF & operator=(const G4StatMF & right);
	G4bool operator==(const G4StatMF & right);
	G4bool operator!=(const G4StatMF & right);

public:

	G4FragmentVector *BreakItUp(const G4Fragment &theNucleus);




private:

	// This finds temperature of breaking channel.
	G4bool FindTemperatureOfBreakingChannel(const G4Fragment & theFragment, 
					  const G4StatMFChannel * aChannel,G4double & Temperature);

	// 
	G4double CalcEnergy(const G4double A, const G4double Z, 
							  const G4StatMFChannel * aChannel,
							  const G4double T);


private:

	G4VStatMFEnsemble * _theEnsemble;



};

#endif
