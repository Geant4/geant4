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
// $Id: G4AlphaGEMProbability.cc,v 1.3 2005/06/04 13:25:25 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4AlphaGEMProbability.hh"

G4AlphaGEMProbability::G4AlphaGEMProbability() :
    G4GEMProbability(4,2,0.0) // A,Z,Gamma
{
    ExcitEnergies.push_back(20.01E+3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(207.0*keV);

    ExcitEnergies.push_back(21.18E+3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(0.73*MeV);

    ExcitEnergies.push_back(22.02E+3*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(1.83*MeV);

    ExcitEnergies.push_back(25.33E+3*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(2.36*MeV);

    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4AlphaGEMProbability::G4AlphaGEMProbability(const G4AlphaGEMProbability &) : G4GEMProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4AlphaGEMProbability::copy_constructor meant to not be accessable");
}




const G4AlphaGEMProbability & G4AlphaGEMProbability::
operator=(const G4AlphaGEMProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4AlphaGEMProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4AlphaGEMProbability::operator==(const G4AlphaGEMProbability &) const
{
    return false;
}

G4bool G4AlphaGEMProbability::operator!=(const G4AlphaGEMProbability &) const
{
    return true;
}


G4double G4AlphaGEMProbability::CCoeficient(const G4double aZ) const
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

