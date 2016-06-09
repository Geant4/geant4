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
// $Id: G4NeutronEvaporationProbability.cc,v 1.4 2006/06/29 20:10:39 gunter Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4NeutronEvaporationProbability.hh"

G4NeutronEvaporationProbability::G4NeutronEvaporationProbability() :
    G4EvaporationProbability(1,0,2) // A,Z,Gamma
{
    std::vector<G4double>::size_type NumExcitedStatesEnergy = 31+1;
    std::vector<G4int>::size_type NumExcitedStatesSpin = 31+1;
    ExcitEnergies.reserve(NumExcitedStatesEnergy);
    ExcitSpins.reserve(NumExcitedStatesSpin);
    ExcitEnergies.insert(ExcitEnergies.begin(),NumExcitedStatesEnergy,0.0);
    ExcitSpins.insert(ExcitSpins.begin(),NumExcitedStatesSpin,0);

	
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


G4NeutronEvaporationProbability::G4NeutronEvaporationProbability(const G4NeutronEvaporationProbability &) : G4EvaporationProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4NeutronEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4NeutronEvaporationProbability & G4NeutronEvaporationProbability::
operator=(const G4NeutronEvaporationProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4NeutronEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4NeutronEvaporationProbability::operator==(const G4NeutronEvaporationProbability &) const
{
    return false;
}

G4bool G4NeutronEvaporationProbability::operator!=(const G4NeutronEvaporationProbability &) const
{
    return true;
}
