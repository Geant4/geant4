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
// $Id: G4AlphaEvaporationProbability.cc,v 1.4 2006/06/29 20:10:19 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4AlphaEvaporationProbability.hh"

G4AlphaEvaporationProbability::G4AlphaEvaporationProbability() :
    G4EvaporationProbability(4,2,4) // A,Z,Gamma
{
    //  const G4int NumExcitedStates = 31+1;
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0);
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


G4AlphaEvaporationProbability::G4AlphaEvaporationProbability(const G4AlphaEvaporationProbability &): G4EvaporationProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4AlphaEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4AlphaEvaporationProbability & G4AlphaEvaporationProbability::
operator=(const G4AlphaEvaporationProbability &) 
{
    throw G4HadronicException(__FILE__, __LINE__, "G4AlphaEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4AlphaEvaporationProbability::operator==(const G4AlphaEvaporationProbability &) const
{
    return false;
}

G4bool G4AlphaEvaporationProbability::operator!=(const G4AlphaEvaporationProbability &) const
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
