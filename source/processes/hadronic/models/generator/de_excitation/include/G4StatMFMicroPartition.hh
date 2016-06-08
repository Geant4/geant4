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
// $Id: G4StatMFMicroPartition.hh,v 1.9 2002/06/18 14:13:01 gcosmo Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#ifndef G4StatMFMicroPartition_h
#define G4StatMFMicroPartition_h 1

#include  "g4std/vector"

#include "globals.hh"
#include "G4StatMFParameters.hh"
#include "G4StatMFChannel.hh"

class G4StatMFMicroPartition {

public:
    // Constructor
    G4StatMFMicroPartition(const G4int A, const G4double Z) :
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
    G4StatMFChannel * ChooseZ(const G4double A0, const G4double Z0, const G4double MeanT);	

    G4double GetProbability(void)
	{ return _Probability; }
	


    void SetPartitionFragment(const G4int anA)
	{ 
	    _thePartition.push_back(anA);
	    CoulombFreeEnergy(anA);
	}
	
    void Normalize(const G4double Normalization)
	{ _Probability /= Normalization; }


    G4double CalcPartitionProbability(const G4double U,
				      const G4double FreeInternalE0,
				      const G4double SCompound);
								
    G4double GetTemperature(void)
	{
	    return _Temperature;
	}
	
    G4double GetEntropy(void)
	{
	    return _Entropy;
	}

private:
    void CoulombFreeEnergy(const G4double anA);

    G4double CalcPartitionTemperature(const G4double U, 
				      const G4double FreeInternalE0);

    G4double GetPartitionEnergy(const G4double T);
	
    G4double GetCoulombEnergy(void);

    G4double GetDegeneracyFactor(const G4int A);
	
    G4double InvLevelDensity(const G4double Af) 
	{
	    // Calculate Inverse Density Level
	    // Epsilon0*(1 + 3 /(Af - 1))
	    if (Af < 1.5) return 0.0;
	    else return G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(Af - 1.0));
	}

private:

    // A and Z of initial nucleus
    G4double theA;
    G4double theZ;

    // Partition probability
    G4double _Probability;
	
    // Partition temperature
    G4double _Temperature;
	
    // Partition entropy
    G4double _Entropy;
	
    // The partition itself
    G4std::vector<G4int> _thePartition;
	
	
    G4std::vector<G4double> _theCoulombFreeEnergy;

};

#endif
