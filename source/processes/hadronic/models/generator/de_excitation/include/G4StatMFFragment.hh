#ifndef G4StatMFFragment_h
#define G4StatMFFragment_h 1

#include "G4StatMFParameters.hh"
#include "G4ThreeVector.hh"
#include "G4Fragment.hh"

class G4StatMFFragment {

public:
	// Constructor
	G4StatMFFragment(const G4int anA, const G4int aZ) :
	theA(anA),theZ(aZ),
	_position(0.0,0.0,0.0),
	_momentum(0.0,0.0,0.0)
	{}


	// Destructor
	virtual ~G4StatMFFragment() {};


private:
	// Default constructor
	G4StatMFFragment(){};
	
	// Copy constructor
	G4StatMFFragment(const G4StatMFFragment & right);

	// operators
	G4StatMFFragment & operator=(const G4StatMFFragment & right);
public:
	G4bool operator==(const G4StatMFFragment & right) const;
	G4bool operator!=(const G4StatMFFragment & right) const;
	
public:

	G4double GetCoulombEnergy(void);
	
	G4double GetEnergy(const G4double T);
	
	G4double GetInvLevelDensity(void);

	G4double GetA(void) const {return theA;}
	
	G4double GetZ(void) const {return theZ;}
	
	void SetPosition(const G4ThreeVector aPosition) {_position = aPosition;}
	
	G4ThreeVector GetPosition(void) {return _position;}
	
	void SetMomentum(const G4ThreeVector aMomentum) {_momentum = aMomentum;}

	G4ThreeVector GetMomentum(void) {return _momentum;}

	G4Fragment * GetFragment(const G4double T);
	
	G4double GetNuclearMass(void)
	{return G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(theZ,theA);}
	

private:

	G4double CalcExcitationEnergy(const G4double T);

private:

	G4double theA;
	
	G4double theZ;
	
	G4ThreeVector _position;
	
	G4ThreeVector _momentum;
};

#endif
