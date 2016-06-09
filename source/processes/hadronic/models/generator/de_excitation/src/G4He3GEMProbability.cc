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
// $Id: G4He3GEMProbability.cc,v 1.1 2002/06/06 18:01:32 larazb Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4He3GEMProbability.hh"

G4He3GEMProbability::G4He3GEMProbability() :
    G4GEMProbability(3,2,1.0/2.0) // A,Z,Gamma
{
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4He3GEMProbability::G4He3GEMProbability(const G4He3GEMProbability &right)
{
    G4Exception("G4He3GEMProbability::copy_constructor meant to not be accessable");
}




const G4He3GEMProbability & G4He3GEMProbability::
operator=(const G4He3GEMProbability &right)
{
    G4Exception("G4He3GEMProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4He3GEMProbability::operator==(const G4He3GEMProbability &right) const
{
    return false;
}

G4bool G4He3GEMProbability::operator!=(const G4He3GEMProbability &right) const
{
    return true;
}


G4double G4He3GEMProbability::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    //	G4double Calpha[5] = { 0.10, 0.10, 0.10, 0.08, 0.06};
    // C for He3 is equal to C for alpha times 4/3
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
    return C*(4.0/3.0);
}

