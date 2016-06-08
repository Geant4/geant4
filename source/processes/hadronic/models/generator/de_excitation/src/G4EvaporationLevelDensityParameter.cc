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
// $Id: G4EvaporationLevelDensityParameter.cc,v 1.4.2.1 2001/06/28 19:13:15 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4EvaporationLevelDensityParameter.hh"


const G4double G4EvaporationLevelDensityParameter::ConstEvapLevelDensityParameter = 0.125*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::alpha = 0.072*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::beta = 0.275*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::gamma = 0.052*(1./MeV);
const G4double G4EvaporationLevelDensityParameter::Bs = 1.0;


G4EvaporationLevelDensityParameter::
G4EvaporationLevelDensityParameter(const G4EvaporationLevelDensityParameter &right)
{
    G4Exception("G4EvaporationLevelDensityParameter::copy_constructor meant to not be accessable");
}


const G4EvaporationLevelDensityParameter & G4EvaporationLevelDensityParameter::
operator=(const G4EvaporationLevelDensityParameter &right)
{
    G4Exception("G4EvaporationLevelDensityParameter::operator= meant to not be accessable");
    return *this;
}


G4bool G4EvaporationLevelDensityParameter::operator==(const G4EvaporationLevelDensityParameter &right) const
{
    return false;
}

G4bool G4EvaporationLevelDensityParameter::operator!=(const G4EvaporationLevelDensityParameter &right) const
{
    return true;
}

G4double G4EvaporationLevelDensityParameter::LevelDensityParameter(const G4int A,const G4int Z,
								   const G4double U) const 
{
    G4int N = A - Z;

    // Asymptotic Level Density Parameter
    G4double AsymptoticLDP = (alpha*G4double(A) + beta*pow(G4double(A),2./3.)*Bs)/MeV;
	
    // Shape of the LDP U dependence
    G4double exponent = -gamma*U;
    G4double f = 1.;
    if (exponent > -300.) f -= exp(exponent);
	
    // Level Density Parameter
    G4double a = AsymptoticLDP*(1. + ShellCorrection(Z,N)*f/U);
	
    return a;
}

