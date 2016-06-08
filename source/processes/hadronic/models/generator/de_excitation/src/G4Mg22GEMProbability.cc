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
// $Id: G4Mg22GEMProbability.cc,v 1.1 2002/06/06 18:02:09 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Mg22GEMProbability.hh"

G4Mg22GEMProbability::G4Mg22GEMProbability() :
  G4GEMProbability(22,12,0.0) // A,Z,Spin
{

    ExcitEnergies.push_back(1246.3*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(2.1*picosecond);

    ExcitEnergies.push_back(3308.2*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(201.0e-3*picosecond);

    ExcitEnergies.push_back(4400.9*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(21.0e-3*picosecond);

    ExcitEnergies.push_back(5006.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(17.0*nanosecond);

    ExcitEnergies.push_back(5037.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(69.0*picosecond);

    ExcitEnergies.push_back(5292.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(44.0e-3*picosecond);

    ExcitEnergies.push_back(5317.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(17.0*nanosecond);

    ExcitEnergies.push_back(5464.0*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(69.0*picosecond);

    ExcitEnergies.push_back(5713.9*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(28.0e-3*picosecond);

}


G4Mg22GEMProbability::G4Mg22GEMProbability(const G4Mg22GEMProbability &right)
{
  G4Exception("G4Mg22GEMProbability::copy_constructor meant to not be accessable");
}




const G4Mg22GEMProbability & G4Mg22GEMProbability::
operator=(const G4Mg22GEMProbability &right)
{
  G4Exception("G4Mg22GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Mg22GEMProbability::operator==(const G4Mg22GEMProbability &right) const
{
  return false;
}

G4bool G4Mg22GEMProbability::operator!=(const G4Mg22GEMProbability &right) const
{
  return true;
}



