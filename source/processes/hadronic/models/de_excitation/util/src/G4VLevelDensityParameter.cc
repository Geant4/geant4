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
// $Id: G4VLevelDensityParameter.cc,v 1.2 2003/11/03 17:53:06 hpw Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4VLevelDensityParameter.hh"


G4VLevelDensityParameter::
G4VLevelDensityParameter(const G4VLevelDensityParameter &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VLevelDensityParameter::copy_constructor meant to not be accessable");
}




const G4VLevelDensityParameter & G4VLevelDensityParameter::
operator=(const G4VLevelDensityParameter &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VLevelDensityParameter::operator= meant to not be accessable");
    return *this;
}


G4bool G4VLevelDensityParameter::
operator==(const G4VLevelDensityParameter &) const
{
    return false;
}

G4bool G4VLevelDensityParameter::
operator!=(const G4VLevelDensityParameter &) const
{
    return true;
}





