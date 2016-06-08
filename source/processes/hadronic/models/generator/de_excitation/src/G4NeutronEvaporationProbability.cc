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
// $Id: G4NeutronEvaporationProbability.cc,v 1.4 2001/08/01 17:05:31 hpw Exp $
// GEANT4 tag $Name: geant4-04-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4NeutronEvaporationProbability.hh"

G4NeutronEvaporationProbability::G4NeutronEvaporationProbability() :
    G4EvaporationProbability(1,0,2) // A,Z,Gamma
{
    G4std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    G4std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0.0);

	
    ExcitEnergies[ 9] = 3.56*MeV;
    ExcitEnergies[10] = 0.48*MeV;
    ExcitEnergies[11] = 0.98*MeV;
    ExcitEnergies[12] = 0.43*MeV;
    ExcitEnergies[15] = 3.37*MeV;
    ExcitEnergies[17] = 0.72*MeV;
    ExcitEnergies[18] = 2.13*MeV;
    ExcitEnergies[19] = 0.95*MeV;
    ExcitEnergies[20] = 2.00*MeV;
    ExcitEnergies[21] = 4.44*MeV;
    ExcitEnergies[22] = 3.09*MeV;
    ExcitEnergies[23] = 6.09*MeV;
    ExcitEnergies[25] = 2.31*MeV;
    ExcitEnergies[26] = 5.28*MeV;
    ExcitEnergies[27] = 0.12*MeV;
    ExcitEnergies[28] = 5.22*MeV;
    ExcitEnergies[29] = 6.10*MeV;
    ExcitEnergies[30] = 0.87*MeV;
    ExcitEnergies[31] = 1.98*MeV;


    ExcitSpins[ 9] = 1;
    ExcitSpins[10] = 2;
    ExcitSpins[11] = 3;
    ExcitSpins[12] = 2;
    ExcitSpins[15] = 5;
    ExcitSpins[17] = 3;
    ExcitSpins[18] = 2;
    ExcitSpins[19] = 5;
    ExcitSpins[20] = 2;
    ExcitSpins[21] = 5;
    ExcitSpins[22] = 2;
    ExcitSpins[23] = 3;
    ExcitSpins[25] = 1;
    ExcitSpins[26] = 8;
    ExcitSpins[27] = 1;
    ExcitSpins[28] = 8;
    ExcitSpins[29] = 8;
    ExcitSpins[30] = 2;
    ExcitSpins[31] = 5;

    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
}


G4NeutronEvaporationProbability::G4NeutronEvaporationProbability(const G4NeutronEvaporationProbability &right)
{
    G4Exception("G4NeutronEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4NeutronEvaporationProbability & G4NeutronEvaporationProbability::
operator=(const G4NeutronEvaporationProbability &right)
{
    G4Exception("G4NeutronEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4NeutronEvaporationProbability::operator==(const G4NeutronEvaporationProbability &right) const
{
    return false;
}

G4bool G4NeutronEvaporationProbability::operator!=(const G4NeutronEvaporationProbability &right) const
{
    return true;
}
