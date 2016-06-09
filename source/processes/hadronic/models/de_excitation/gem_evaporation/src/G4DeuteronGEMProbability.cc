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
// $Id: G4DeuteronGEMProbability.cc,v 1.2 2003/11/03 17:53:03 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4DeuteronGEMProbability.hh"

G4DeuteronGEMProbability::G4DeuteronGEMProbability() :
    G4GEMProbability(2,1,1.0) // A,Z,Gamma
{
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4DeuteronGEMProbability::G4DeuteronGEMProbability(const G4DeuteronGEMProbability &) : G4GEMProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4DeuteronGEMProbability::copy_constructor meant to not be accessable");
}




const G4DeuteronGEMProbability & G4DeuteronGEMProbability::
operator=(const G4DeuteronGEMProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4DeuteronGEMProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4DeuteronGEMProbability::operator==(const G4DeuteronGEMProbability &) const
{
    return false;
}

G4bool G4DeuteronGEMProbability::operator!=(const G4DeuteronGEMProbability &) const
{
    return true;
}


G4double G4DeuteronGEMProbability::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
    // C for deuteron is equal to C for protons divided by 2
    G4double C = 0.0;
	
    if (aZ >= 70) {
        C = 0.10;
    } else {
        C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
    }
	
    return C/2.0;
	
}
