//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StatMFMacroTemperature.hh,v 1.1 2003/08/26 18:47:21 lara Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroTemperature_h
#define G4StatMFMacroTemperature_h 1

#include "G4StatMFParameters.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4StatMFMacroChemicalPotential.hh"
#include "G4Solver.hh"



class G4StatMFMacroTemperature {

public:

    G4StatMFMacroTemperature(const G4double anA, const G4double aZ, 
			     const G4double ExEnergy, const G4double FreeE0, 
			     const G4double kappa, 
			     std::vector<G4VStatMFMacroCluster*> * ClusterVector) :
	theA(anA),
	theZ(aZ),
	_ExEnergy(ExEnergy),
	_FreeInternalE0(FreeE0),
	_Kappa(kappa),
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
    G4StatMFMacroTemperature(const G4StatMFMacroTemperature &) {};


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
	
    std::vector<G4VStatMFMacroCluster*> * _theClusters; 


};
#endif
