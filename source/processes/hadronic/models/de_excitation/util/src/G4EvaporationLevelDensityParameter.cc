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
// $Id: G4EvaporationLevelDensityParameter.cc,v 1.4 2005/06/04 13:29:20 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4EvaporationLevelDensityParameter.hh"
#include "G4HadronicException.hh"

// Those values are from table 3 in 
// A.S. Iljinov et al. Nucl Phys A543 (1992) 517-557
// Table 3. alpha, beta and gamma for Cameron Shell corrections
// whithout collective effects. f-factor = 2.31.

const G4double G4EvaporationLevelDensityParameter::ConstEvapLevelDensityParameter = 0.125*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::alpha = 0.072*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::beta = 0.257*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::gamma = 0.059*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::Bs = 1.0;


G4EvaporationLevelDensityParameter::
G4EvaporationLevelDensityParameter(const G4EvaporationLevelDensityParameter &) : G4VLevelDensityParameter()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationLevelDensityParameter::copy_constructor meant to not be accessable");
}


const G4EvaporationLevelDensityParameter & G4EvaporationLevelDensityParameter::
operator=(const G4EvaporationLevelDensityParameter &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4EvaporationLevelDensityParameter::operator= meant to not be accessable");
    return *this;
}


G4bool G4EvaporationLevelDensityParameter::operator==(const G4EvaporationLevelDensityParameter &) const
{
    return false;
}

G4bool G4EvaporationLevelDensityParameter::operator!=(const G4EvaporationLevelDensityParameter &) const
{
    return true;
}

G4double G4EvaporationLevelDensityParameter::LevelDensityParameter(const G4int A,const G4int Z,
								   const G4double U) const 
{
    G4int N = A - Z;

    // Asymptotic Level Density Parameter
    G4double AsymptoticLDP = (alpha*static_cast<G4double>(A) + beta*std::pow(static_cast<G4double>(A),2./3.)*Bs)/MeV;
	
    // Shape of the LDP U dependence
    G4double exponent = -gamma*U;
    G4double f = 1.;
    if (exponent > -300.) f -= std::exp(exponent);
	
    // Level Density Parameter
    G4double a = AsymptoticLDP*(1. + ShellCorrection(Z,N)*f/U);
	
    return a;
}

