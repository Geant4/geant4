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
// $Id: G4DeuteronEvaporationProbability.cc,v 1.4 2006/06/29 20:10:23 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4DeuteronEvaporationProbability.hh"

G4DeuteronEvaporationProbability::G4DeuteronEvaporationProbability() :
    G4EvaporationProbability(2,1,6) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0);


    ExcitEnergies[15] = 6.18*MeV;
    ExcitEnergies[17] = 2.15*MeV;
    ExcitEnergies[18] = 5.02*MeV;
    ExcitEnergies[19] = 2.65*MeV;
    ExcitEnergies[20] = 4.80*MeV;
    ExcitEnergies[22] = 3.85*MeV;
    ExcitEnergies[23] = 6.96*MeV;
    ExcitEnergies[25] = 4.92*MeV;
    ExcitEnergies[26] = 7.22*MeV;
    ExcitEnergies[27] = 0.40*MeV;
    ExcitEnergies[28] = 6.83*MeV;
    ExcitEnergies[29] = 7.12*MeV;
    ExcitEnergies[30] = 3.84*MeV;
    ExcitEnergies[31] = 3.92*MeV;

    ExcitSpins[15] = 1;
    ExcitSpins[17] = 3;
    ExcitSpins[18] = 4;
    ExcitSpins[19] = 4;
    ExcitSpins[20] = 4;
    ExcitSpins[22] = 6;
    ExcitSpins[23] = 6;
    ExcitSpins[25] = 1;
    ExcitSpins[26] = 10;
    ExcitSpins[27] = 3;
    ExcitSpins[28] = 10;
    ExcitSpins[29] = 3;
    ExcitSpins[30] = 6;
    ExcitSpins[31] = 5;
	
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);	
	
}
G4DeuteronEvaporationProbability::G4DeuteronEvaporationProbability(const G4DeuteronEvaporationProbability &) : G4EvaporationProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4DeuteronEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4DeuteronEvaporationProbability & G4DeuteronEvaporationProbability::
operator=(const G4DeuteronEvaporationProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4DeuteronEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4DeuteronEvaporationProbability::operator==(const G4DeuteronEvaporationProbability &) const
{
    return false;
}

G4bool G4DeuteronEvaporationProbability::operator!=(const G4DeuteronEvaporationProbability &) const
{
    return true;
}

G4double G4DeuteronEvaporationProbability::CCoeficient(const G4double aZ) const
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
