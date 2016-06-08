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
// $Id: G4FissionLevelDensityParameter.cc,v 1.6 2001/08/01 17:05:31 hpw Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4FissionLevelDensityParameter.hh"


G4FissionLevelDensityParameter::
G4FissionLevelDensityParameter(const G4FissionLevelDensityParameter &right)
{
    G4Exception("G4FissionLevelDensityParameter::copy_constructor meant to not be accessable");
}


const G4FissionLevelDensityParameter & G4FissionLevelDensityParameter::
operator=(const G4FissionLevelDensityParameter &right)
{
    G4Exception("G4FissionLevelDensityParameter::operator= meant to not be accessable");
    return *this;
}


G4bool G4FissionLevelDensityParameter::
operator==(const G4FissionLevelDensityParameter &right) const
{
    return false;
}

G4bool G4FissionLevelDensityParameter::
operator!=(const G4FissionLevelDensityParameter &right) const
{
    return true;
}


G4double G4FissionLevelDensityParameter::
LevelDensityParameter(const G4int A,const G4int Z,const G4double U) const 
{
    G4double EvapLDP = theEvaporationLevelDensityParameter.LevelDensityParameter(A,Z,U);

    if (Z >= 89) return 1.04*EvapLDP;
    else if (Z >= 85) return (1.04*(1./MeV) + 0.01*(89-Z))*EvapLDP;
    else return 1.08*EvapLDP;
}
