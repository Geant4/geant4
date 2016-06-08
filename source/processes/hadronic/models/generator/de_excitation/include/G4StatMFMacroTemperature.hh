#ifndef G4StatMFMacroTemperature_h
#define G4StatMFMacroTemperature_h 1

#include "g4rw/tpordvec.h"

#include "G4StatMFParameters.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4StatMFMacroChemicalPotential.hh"
#include "G4Solver.hh"



class G4StatMFMacroTemperature {

public:

	G4StatMFMacroTemperature(const G4double anA, const G4double aZ, 
									 const G4double ExEnergy, const G4double FreeE0, 
									 const G4double kappa, 
									 G4RWTPtrOrderedVector<G4VStatMFMacroCluster> * ClusterVector) :
		theA(anA),
		theZ(aZ),
		_ExEnergy(ExEnergy),
		_Kappa(kappa),
		_FreeInternalE0(FreeE0),
		_MeanMultiplicity(0.0),
		_MeanTemperature(0.0),
		_ChemPotentialMu(0.0),
		_ChemPotentialNu(0.0),
		_theClusters(ClusterVector) 
		{};
	
	~G4StatMFMacroTemperature() {};
   
	G4double operator()(const G4double T)
	{ return (_ExEnergy - this->FragsExcitEnergy(T))/_ExEnergy; }	

private:
	// Default constructor
	G4StatMFMacroTemperature() {};

	// copy constructor
	G4StatMFMacroTemperature(const G4StatMFMacroTemperature &right) {};


	// operators
	G4StatMFMacroTemperature & operator=(const G4StatMFMacroTemperature & right);
	G4bool operator==(const G4StatMFMacroTemperature & right) const;
	G4bool operator!=(const G4StatMFMacroTemperature & right) const;

public:

	G4double GetMeanMultiplicity(void) const {return _MeanMultiplicity;}
	
	G4double GetChemicalPotentialMu(void) const {return _ChemPotentialMu;}

	G4double GetChemicalPotentialNu(void) const {return _ChemPotentialNu;}

	G4double GetTemperature(void) const {return _MeanTemperature;}

	G4double GetEntropy(void) const {return _MeanEntropy;}

	G4double CalcTemperature(void);

private:
	
	G4double FragsExcitEnergy(const G4double T);

	void CalcChemicalPotentialNu(const G4double T);

private:

	G4double theA;

	G4double theZ;

	G4double _ExEnergy;
	
	G4double _FreeInternalE0;

	G4double _Kappa;

	G4double _MeanMultiplicity;

	G4double _MeanTemperature;
	
	G4double _ChemPotentialMu;
	
	G4double _ChemPotentialNu;
	
	G4double _MeanEntropy;
	
	G4RWTPtrOrderedVector<G4VStatMFMacroCluster> * _theClusters; 


};
#endif
