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
// $Id: G4ProtonGEMProbability.cc,v 1.1 2002/06/06 18:03:33 larazb Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Sept 2001)
//


#include "G4ProtonGEMProbability.hh"

G4ProtonGEMProbability::G4ProtonGEMProbability() :
    G4GEMProbability(1,1,1.0/2.0) // A,Z,Gamma
{
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4ProtonGEMProbability::G4ProtonGEMProbability(const G4ProtonGEMProbability &right)
{
    G4Exception("G4ProtonGEMProbability::copy_constructor meant to not be accessable");
}




const G4ProtonGEMProbability & G4ProtonGEMProbability::
operator=(const G4ProtonGEMProbability &right)
{
    G4Exception("G4ProtonGEMProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4ProtonGEMProbability::operator==(const G4ProtonGEMProbability &right) const
{
    return false;
}

G4bool G4ProtonGEMProbability::operator!=(const G4ProtonGEMProbability &right) const
{
    return true;
}


G4double G4ProtonGEMProbability::CCoeficient(const G4double aZ) const
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
