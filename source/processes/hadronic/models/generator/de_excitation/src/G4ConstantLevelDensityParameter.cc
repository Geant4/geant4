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
// Hadronic Process: Nuclear De-excitations
//    Constant level density parameter (for photon evaporation)
//
// by C. Dallapiccola (Nov 1998)
//


#include "G4ConstantLevelDensityParameter.hh"

G4ConstantLevelDensityParameter::
G4ConstantLevelDensityParameter(const G4ConstantLevelDensityParameter& ) :
 G4VLevelDensityParameter(),  EvapLevelDensityParameter(0.125*(1./MeV))
{
  G4Exception("G4ConstantLevelDensityParameter::copy_constructor meant to not be accessable");
}


const G4ConstantLevelDensityParameter & G4ConstantLevelDensityParameter::
operator=(const G4ConstantLevelDensityParameter &)
{
  G4Exception("G4ConstantLevelDensityParameter::operator= meant to not be accessable");
  return *this;
}


G4bool G4ConstantLevelDensityParameter::operator==(const G4ConstantLevelDensityParameter &) const
{
  return false;
}

G4bool G4ConstantLevelDensityParameter::operator!=(const G4ConstantLevelDensityParameter &) const
{
  return true;
}





