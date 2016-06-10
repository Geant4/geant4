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
// $Id: G4StatMFMicroPartition.hh 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMicroPartition_h
#define G4StatMFMicroPartition_h 1

#include  <vector>

#include "globals.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"

class G4StatMFMicroPartition {

public:
    // Constructor
    G4StatMFMicroPartition(G4int A, G4int Z) :
	theA(A), theZ(Z), _Probability(0.0), _Temperature(0.0), 
	_Entropy(0.0) {};


    // Destructor
    ~G4StatMFMicroPartition() {};


private:
    // Default constructor
    G4StatMFMicroPartition() {};
	
    // Copy constructor
    G4StatMFMicroPartition(const G4StatMFMicroPartition & right);

    // operators
    G4StatMFMicroPartition & operator=(const G4StatMFMicroPartition & right);
public:
    G4bool operator==(const G4StatMFMicroPartition & right) const;
    G4bool operator!=(const G4StatMFMicroPartition & right) const;

public:

    // Gives fragments charges
    G4StatMFChannel * ChooseZ(G4int A0, G4int Z0, G4double MeanT);	

    G4double GetProbability(void)
	{ return _Probability; }
	
    void SetPartitionFragment(G4int anA)
	{ 
	    _thePartition.push_back(anA);
	    CoulombFreeEnergy(anA);
	}
	
    void Normalize(G4double Normalization)
	{ _Probability /= Normalization; }

    G4double CalcPartitionProbability(G4double U,
				      G4double FreeInternalE0,
				      G4double SCompound);
								
    G4double GetTemperature(void)
	{
	    return _Temperature;
	}
	
    G4double GetEntropy(void)
	{
	    return _Entropy;
	}

private:

    void CoulombFreeEnergy(G4int anA);

    G4double CalcPartitionTemperature(G4double U, 
				      G4double FreeInternalE0);

    G4double GetPartitionEnergy(G4double T);
	
    G4double GetCoulombEnergy(void);

    G4double GetDegeneracyFactor(G4int A);
	
    G4double InvLevelDensity(G4double Af) 
	{
	    // Calculate Inverse Density Level
	    // Epsilon0*(1 + 3 /(Af - 1))
	    if (Af < 1.5) return 0.0;
	    else return G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(Af - 1.0));
	}

private:

    // A and Z of initial nucleus
    G4int theA;
    G4int theZ;

    // Partition probability
    G4double _Probability;
	
    // Partition temperature
    G4double _Temperature;
	
    // Partition entropy
    G4double _Entropy;
	
    // The partition itself
    std::vector<G4int> _thePartition;
	
    std::vector<G4double> _theCoulombFreeEnergy;

};

#endif
