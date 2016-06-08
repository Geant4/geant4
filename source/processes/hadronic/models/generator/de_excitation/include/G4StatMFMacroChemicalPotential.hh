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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4StatMFMacroChemicalPotential.hh,v 1.7 2002/06/18 14:13:01 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroChemicalPotential_h
#define G4StatMFMacroChemicalPotential_h 1

#include "g4std/vector"

#include "G4StatMFParameters.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4StatMFMacroMultiplicity.hh"
#include "G4Solver.hh"



class G4StatMFMacroChemicalPotential {

public:

    G4StatMFMacroChemicalPotential(const G4double anA, const G4double aZ,
				   const G4double kappa, 
				   const G4double temp, 
				   G4std::vector<G4VStatMFMacroCluster*> * ClusterVector) :
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
	
    G4std::vector<G4VStatMFMacroCluster*> * _theClusters; 


};
#endif
