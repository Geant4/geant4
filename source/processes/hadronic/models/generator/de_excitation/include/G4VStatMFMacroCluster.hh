#ifndef G4VStatMFMacroCluster_h
#define G4VStatMFMacroCluster_h 1

#include "G4StatMFParameters.hh"

class G4VStatMFMacroCluster {

public:
	// Constructor
	G4VStatMFMacroCluster(const G4int Size) : 
		theA(Size),
	_InvLevelDensity(0.0),
	_Entropy(0.0),
	theZARatio(0.0),
	_MeanMultiplicity(0.0),
	_Energy(0.0)
	{
		if (theA <= 0) G4Exception(
		"G4VStatMFMacroCluster::Constructor: Cluster's size must be >= 1");
	 	_InvLevelDensity = CalcInvLevelDensity();
	 }


	// Destructor
	virtual ~G4VStatMFMacroCluster() {};


private:

	// Default constructor
	G4VStatMFMacroCluster() {};

	// Copy constructor
	G4VStatMFMacroCluster(const G4VStatMFMacroCluster & right);

	// operators
	G4VStatMFMacroCluster & operator=(const G4VStatMFMacroCluster & right);

public:
	G4bool operator==(const G4VStatMFMacroCluster & right) const;
	G4bool operator!=(const G4VStatMFMacroCluster & right) const;

private:
	G4double CalcInvLevelDensity(void);
	
public:

	virtual G4double CalcMeanMultiplicity(const G4double FreeVol, const G4double mu, 
													const G4double nu, const G4double T) = 0;
						
	virtual G4double CalcZARatio(const G4double nu) = 0;
	
	G4double GetMeanMultiplicity(void) const { return _MeanMultiplicity; }

	virtual G4double CalcEnergy(const G4double T) = 0;

	virtual G4double CalcEntropy(const G4double T, const G4double FreeVol) = 0;

protected:
	// Number of nucleons in the cluster
	G4int theA;
	
	// Inverse level density
	G4double _InvLevelDensity;
	
	// Entropy
	G4double _Entropy;
	
	// Z/A ratio
	G4double theZARatio;
	
	// Mean Multiplicity
	G4double _MeanMultiplicity;
	
	// Energy
	G4double _Energy;

// *************************************************************************
public:
	
	const G4double GetInvLevelDensity(void) const
	{ return _InvLevelDensity; }
	

	 
	void SetZARatio(const G4double value)
	{ theZARatio = value; }
	
	const G4double GetZARatio(void) const
	{ return theZARatio; }
	
	
	void SetSize(const G4double value)
	{ 
		if (value <= 0.0) G4Exception("G4VStatMFMacroCluster::SetSize: Cluster's size must be >= 1");
		theA = G4int(value); 
		_InvLevelDensity = CalcInvLevelDensity();
	}
	
	const G4double GetSize(void) const
	{ return theA; }





};

#endif
