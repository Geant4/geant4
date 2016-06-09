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
// $Id: G4FissionProbability.cc,v 1.8 2007/02/12 09:39:58 ahoward Exp $
// GEANT4 tag $Name: geant4-09-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4FissionProbability.hh"
#include "G4PairingCorrection.hh"



G4FissionProbability::G4FissionProbability(const G4FissionProbability &) : G4VEmissionProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4FissionProbability::copy_constructor meant to not be accessable");
}




const G4FissionProbability & G4FissionProbability::operator=(const G4FissionProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4FissionProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4FissionProbability::operator==(const G4FissionProbability &) const
{
    return false;
}

G4bool G4FissionProbability::operator!=(const G4FissionProbability &) const
{
    return true;
}


G4double G4FissionProbability::EmissionProbability(const G4Fragment & fragment, const G4double MaximalKineticEnergy)
    // Compute integrated probability of fission channel
{
    if (MaximalKineticEnergy <= 0.0) return 0.0;
    G4double A = fragment.GetA();
    G4double Z = fragment.GetZ();
    G4double U = fragment.GetExcitationEnergy();
  
    G4double Ucompound = U - G4PairingCorrection::GetInstance()->GetPairingCorrection(static_cast<G4int>(A),
										      static_cast<G4int>(Z));
    G4double Ufission = U - G4PairingCorrection::GetInstance()->GetFissionPairingCorrection(static_cast<G4int>(A),
											    static_cast<G4int>(Z));
  
    G4double SystemEntropy = 2.0*std::sqrt(theEvapLDP.LevelDensityParameter(static_cast<G4int>(A),
								       static_cast<G4int>(Z),
								       Ucompound)*Ucompound);
	
    G4double afission = theFissLDP.LevelDensityParameter(static_cast<G4int>(A),
							 static_cast<G4int>(Z),
							 Ufission);

    G4double Cf = 2.0*std::sqrt(afission*MaximalKineticEnergy);

    //    G4double Q1 = 1.0 + (Cf - 1.0)*std::exp(Cf);
    //    G4double Q2 = 4.0*pi*afission*std::exp(SystemEntropy);
	
    //    G4double probability = Q1/Q2;
 
    
    G4double Exp1 = 0.0;
    if (SystemEntropy <= 160.0) Exp1 = std::exp(-SystemEntropy);
    // @@@@@@@@@@@@@@@@@ hpw changed max to min - cannot notify vicente now since cern mail gave up on me...
    G4double Exp2 = std::exp( std::min(700.0,Cf-SystemEntropy) ); 

    //AH fix from Vincente:    G4double probability = (Exp1 + (1.0-Cf)*Exp2) / 4.0*pi*afission;
    G4double probability = (Exp1 + (Cf-1.0)*Exp2) / 4.0*pi*afission;
    
    return probability;
}

