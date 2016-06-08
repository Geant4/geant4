#ifndef G4StatMFMacroMultiplicity_h
#define G4StatMFMacroMultiplicity_h 1

#include "g4rw/tpordvec.h"

#include "G4StatMFParameters.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4Solver.hh"



class G4StatMFMacroMultiplicity {

public:

	G4StatMFMacroMultiplicity(const G4double anA, 
									  const G4double kappa, 
									  const G4double temp, 
									  const G4double nu,
									  G4RWTPtrOrderedVector<G4VStatMFMacroCluster> * ClusterVector) :
		theA(anA),
		_Kappa(kappa),
		_MeanMultiplicity(0.0),
		_MeanTemperature(temp),
		_ChemPotentialMu(0.0),
		_ChemPotentialNu(nu),
		_theClusters(ClusterVector) 
		{};
	
	~G4StatMFMacroMultiplicity() {};
   
	G4double operator()(const G4double mu)
	{ return (theA - this->CalcMeanA(mu))/theA; }	

private:
	// Default constructor
	G4StatMFMacroMultiplicity() {};

	// copy constructor
	G4StatMFMacroMultiplicity(const G4StatMFMacroMultiplicity &right) {};


	// operators
	G4StatMFMacroMultiplicity & operator=(const G4StatMFMacroMultiplicity & right);
	G4bool operator==(const G4StatMFMacroMultiplicity & right) const;
	G4bool operator!=(const G4StatMFMacroMultiplicity & right) const;

public:

	G4double GetMeanMultiplicity(void) const {return _MeanMultiplicity;}
	
	G4double GetChemicalPotentialMu(void) const {return _ChemPotentialMu;}

	G4double CalcChemicalPotentialMu(void);

private:
	
	G4double CalcMeanA(const G4double mu);

private:

	G4double theA;

	G4double _Kappa;

	G4double _MeanMultiplicity;

	G4double _MeanTemperature;
	
	G4double _ChemPotentialMu;
	
	G4double _ChemPotentialNu;
	
	G4RWTPtrOrderedVector<G4VStatMFMacroCluster> * _theClusters; 


};
#endif
