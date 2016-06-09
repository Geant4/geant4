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
// $Id: G4TritonGEMProbability.cc,v 1.4 2006/06/29 20:23:43 gunter Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4TritonGEMProbability.hh"

G4TritonGEMProbability::G4TritonGEMProbability() :
    G4GEMProbability(3,1,1.0/2.0) // A,Z,Gamma
{
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4TritonGEMProbability::G4TritonGEMProbability(const G4TritonGEMProbability &) : G4GEMProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonGEMProbability::copy_constructor meant to not be accessable");
}




const G4TritonGEMProbability & G4TritonGEMProbability::
operator=(const G4TritonGEMProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4TritonGEMProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4TritonGEMProbability::operator==(const G4TritonGEMProbability &) const
{
    return false;
}

G4bool G4TritonGEMProbability::operator!=(const G4TritonGEMProbability &) const
{
    return true;
}


G4double G4TritonGEMProbability::CCoeficient(const G4double aZ) const
{
    // Data comes from 
    // Dostrovsky, Fraenkel and Friedlander
    // Physical Review, vol 116, num. 3 1959
    // 
    // const G4int size = 5;
    // G4double Zlist[5] = { 10.0, 20.0, 30.0, 50.0, 70.0};
    // G4double Cp[5] = { 0.50, 0.28, 0.20, 0.15, 0.10};
    // C for triton is equal to C for protons divided by 3
    G4double C = 0.0;
	
    if (aZ >= 70) {
	C = 0.10;
    } else {
	C = ((((0.15417e-06*aZ) - 0.29875e-04)*aZ + 0.21071e-02)*aZ - 0.66612e-01)*aZ + 0.98375;
    }
	
    return C/3.0;
	
}
