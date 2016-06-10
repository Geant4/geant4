//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4StatMFMacroChemicalPotential.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMacroChemicalPotential_h
#define G4StatMFMacroChemicalPotential_h 1

#include <vector>

#include "G4StatMFParameters.hh"
#include "G4VStatMFMacroCluster.hh"
#include "G4StatMFMacroMultiplicity.hh"
#include "G4Solver.hh"



class G4StatMFMacroChemicalPotential {

public:

    G4StatMFMacroChemicalPotential(const G4double anA, const G4double aZ,
				   const G4double kappa, 
				   const G4double temp, 
				   std::vector<G4VStatMFMacroCluster*> * ClusterVector) :
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
    G4StatMFMacroChemicalPotential(const G4StatMFMacroChemicalPotential &) {};


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
	
    std::vector<G4VStatMFMacroCluster*> * _theClusters; 


};
#endif
