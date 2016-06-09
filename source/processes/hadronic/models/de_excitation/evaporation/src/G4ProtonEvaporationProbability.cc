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
// $Id: G4ProtonEvaporationProbability.cc,v 1.2 2003/11/03 17:53:02 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4ProtonEvaporationProbability.hh"

G4ProtonEvaporationProbability::G4ProtonEvaporationProbability() :
    G4EvaporationProbability(1,1,2) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0);



    ExcitEnergies[15] = 5.96*MeV;
    ExcitEnergies[17] = 1.74*MeV;
    ExcitEnergies[18] = 4.44*MeV;
    ExcitEnergies[19] = 1.67*MeV;
    ExcitEnergies[20] = 4.32*MeV;
    ExcitEnergies[22] = 3.68*MeV;
    ExcitEnergies[23] = 6.69*MeV;
    ExcitEnergies[25] = 3.95*MeV;
    ExcitEnergies[26] = 6.32*MeV;
    ExcitEnergies[27] = 0.30*MeV;
    ExcitEnergies[28] = 6.18*MeV;
    ExcitEnergies[29] = 6.92*MeV;
    ExcitEnergies[30] = 3.06*MeV;
    ExcitEnergies[31] = 3.57*MeV;


    ExcitSpins[15] = 8;
    ExcitSpins[17] = 1;
    ExcitSpins[18] = 6;
    ExcitSpins[19] = 5;
    ExcitSpins[20] = 6;
    ExcitSpins[22] = 4;
    ExcitSpins[23] = 8;
    ExcitSpins[25] = 3;
    ExcitSpins[26] = 4;
    ExcitSpins[27] = 7;
    ExcitSpins[28] = 4;
    ExcitSpins[29] = 5;
    ExcitSpins[30] = 2;
    ExcitSpins[31] = 10;
	
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);	
	
}

G4ProtonEvaporationProbability::G4ProtonEvaporationProbability(const G4ProtonEvaporationProbability &) : G4EvaporationProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4ProtonEvaporationProbability::copy_constructor meant to not be accessable");
}

const G4ProtonEvaporationProbability & G4ProtonEvaporationProbability::
operator=(const G4ProtonEvaporationProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4ProtonEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4ProtonEvaporationProbability::operator==(const G4ProtonEvaporationProbability &) const
{
    return false;
}

G4bool G4ProtonEvaporationProbability::operator!=(const G4ProtonEvaporationProbability &) const
{
    return true;
}

G4double G4ProtonEvaporationProbability::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
    G4double C = 0.0;
	
    if (aZ >= 70) {
	C = 0.10;
    } else {
	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
    }
	
    return C;
	
}
