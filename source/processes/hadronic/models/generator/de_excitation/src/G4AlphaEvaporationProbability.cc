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
// $Id: G4AlphaEvaporationProbability.cc,v 1.4 2001/08/01 17:05:10 hpw Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4AlphaEvaporationProbability.hh"

G4AlphaEvaporationProbability::G4AlphaEvaporationProbability() :
    G4EvaporationProbability(4,2,4) // A,Z,Gamma
{
    //  const G4int NumExcitedStates = 31+1;
    G4std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    G4std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0.0);
//    for (G4int i = 0; i < NumExcitedStates; i++) {
//      ExcitEnergies(i) = 0.0;
//      ExcitSpins(i) = 0;
//    }

    ExcitEnergies[18] = 7.98*MeV;
    ExcitEnergies[20] = 6.90*MeV;
    ExcitEnergies[25] = 5.83*MeV;
    ExcitEnergies[26] = 8.57*MeV;
    ExcitEnergies[31] = 5.33*MeV;

    ExcitSpins[18] = 4;
    ExcitSpins[20] = 6;
    ExcitSpins[25] = 7;
    ExcitSpins[26] = 4;
    ExcitSpins[31] = 13;
	
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);		
}


G4AlphaEvaporationProbability::G4AlphaEvaporationProbability(const G4AlphaEvaporationProbability &right)
{
    G4Exception("G4AlphaEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4AlphaEvaporationProbability & G4AlphaEvaporationProbability::
operator=(const G4AlphaEvaporationProbability &right)
{
    G4Exception("G4AlphaEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4AlphaEvaporationProbability::operator==(const G4AlphaEvaporationProbability &right) const
{
    return false;
}

G4bool G4AlphaEvaporationProbability::operator!=(const G4AlphaEvaporationProbability &right) const
{
    return true;
}

G4double G4AlphaEvaporationProbability::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    //	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
    G4double C = 0.0;
	
	
    if (aZ <= 30) {
	C = 0.10;
    } else if (aZ <= 50) {
	C = 0.1 + -((aZ-50.)/20.)*0.02;
    } else if (aZ < 70) {
	C = 0.08 + -((aZ-70.)/20.)*0.02;
    } else {
	C = 0.06;
    }
    return C;
}
