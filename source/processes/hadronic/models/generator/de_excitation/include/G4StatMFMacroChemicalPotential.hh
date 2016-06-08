#ifndef G4StatMFMacroChemicalPotential_h
#define G4StatMFMacroChemicalPotential_h 1

#include "g4rw/tpordvec.h"

#include "G4StatMFParameters.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4StatMFMacroMultiplicity.hh"
#include "G4Solver.hh"



class G4StatMFMacroChemicalPotential {

public:

	G4StatMFMacroChemicalPotential(const G4double anA, const G4double aZ,
									  const G4double kappa, 
									  const G4double temp, 
									  G4RWTPtrOrderedVector<G4VStatMFMacroCluster> * ClusterVector) :
		theA(anA),
		theZ(aZ),
		_Kappa(kappa),
		_MeanMultiplicity(0.0),
		_MeanTemperature(temp),
		_ChemPotentialMu(0.0),
		_ChemPotentialNu(0.0),
		_theClusters(ClusterVector) 
		{};
	
	~G4StatMFMacroChemicalPotential() {};
   
	G4double operator()(const G4double nu)
	{ return (theZ - this->CalcMeanZ(nu))/theZ; }	

private:
	// Default constructor
	G4StatMFMacroChemicalPotential() {};

	// copy constructor
	G4StatMFMacroChemicalPotential(const G4StatMFMacroChemicalPotential &right) {};


	// operators
	G4StatMFMacroChemicalPotential & operator=(const G4StatMFMacroChemicalPotential & right);
	G4bool operator==(const G4StatMFMacroChemicalPotential & right) const;
	G4bool operator!=(const G4StatMFMacroChemicalPotential & right) const;

public:

	G4double GetMeanMultiplicity(void) const {return _MeanMultiplicity;}
	
	G4double GetChemicalPotentialMu(void) const {return _ChemPotentialMu;}

	G4double GetChemicalPotentialNu(void) const {return _ChemPotentialNu;}

	G4double CalcChemicalPotentialNu(void);

private:
	
	G4double CalcMeanZ(const G4double nu);

	void CalcChemicalPotentialMu(const G4double nu);

private:

	G4double theA;

	G4double theZ;

	G4double _Kappa;

	G4double _MeanMultiplicity;

	G4double _MeanTemperature;
	
	G4double _ChemPotentialMu;
	
	G4double _ChemPotentialNu;
	
	G4RWTPtrOrderedVector<G4VStatMFMacroCluster> * _theClusters; 


};
#endif
