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
// $Id: G4StatMFMacroMultiplicity.hh,v 1.10 2003/06/16 17:06:04 gunter Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroMultiplicity_h
#define G4StatMFMacroMultiplicity_h 1

#include <vector>

#include "G4StatMFParameters.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4Solver.hh"



class G4StatMFMacroMultiplicity {

public:

    G4StatMFMacroMultiplicity(const G4double anA, 
			      const G4double kappa, 
			      const G4double temp, 
			      const G4double nu,
			      std::vector<G4VStatMFMacroCluster*> * ClusterVector) :
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
    G4StatMFMacroMultiplicity(const G4StatMFMacroMultiplicity &) {};


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
	
    std::vector<G4VStatMFMacroCluster*> * _theClusters; 


};
#endif
