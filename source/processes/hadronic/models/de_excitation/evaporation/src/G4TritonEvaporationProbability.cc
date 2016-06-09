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
// $Id: G4TritonEvaporationProbability.cc,v 1.2 2003/11/03 17:53:02 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4TritonEvaporationProbability.hh"

G4TritonEvaporationProbability::G4TritonEvaporationProbability() :
    G4EvaporationProbability(3,1,6) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0);

	
    ExcitEnergies[15] = 6.26*MeV;
    ExcitEnergies[17] = 3.59*MeV;
    ExcitEnergies[18] = 6.76*MeV;
    ExcitEnergies[20] = 6.34*MeV;
    ExcitEnergies[23] = 7.34*MeV;
    ExcitEnergies[25] = 5.11*MeV;
    ExcitEnergies[26] = 7.57*MeV;
    ExcitEnergies[28] = 7.28*MeV;
    ExcitEnergies[31] = 4.46*MeV;

    ExcitSpins[15] = 5;
    ExcitSpins[17] = 5;
    ExcitSpins[18] = 10;
    ExcitSpins[20] = 2;
    ExcitSpins[23] = 5;
    ExcitSpins[25] = 5;
    ExcitSpins[26] = 8;
    ExcitSpins[28] = 8;
    ExcitSpins[31] = 3;

    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);	
}

G4TritonEvaporationProbability::G4TritonEvaporationProbability(const G4TritonEvaporationProbability &) : G4EvaporationProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4TritonEvaporationProbability & G4TritonEvaporationProbability::
operator=(const G4TritonEvaporationProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4TritonEvaporationProbability::operator==(const G4TritonEvaporationProbability &) const
{
    return false;
}

G4bool G4TritonEvaporationProbability::operator!=(const G4TritonEvaporationProbability &) const
{
    return true;
}

G4double G4TritonEvaporationProbability::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
    // C for triton is equal to C for protons divided by 3
    G4double C = 0.0;
	
    if (aZ >= 70) {
	C = 0.10;
    } else {
	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
    }
	
    return C/3.0;
	
}
