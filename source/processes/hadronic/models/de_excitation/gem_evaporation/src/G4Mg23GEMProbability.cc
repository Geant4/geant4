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
// $Id: G4Mg23GEMProbability.cc,v 1.2 2003/11/03 17:53:04 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Mg23GEMProbability.hh"

G4Mg23GEMProbability::G4Mg23GEMProbability() :
  G4GEMProbability(23,12,3.0/2.0) // A,Z,Spin
{

    ExcitEnergies.push_back(450.70*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(1.25*picosecond);

    ExcitEnergies.push_back(2051.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(55.0e-3*picosecond);

    ExcitEnergies.push_back(2359.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(575.0e-3*picosecond);

    ExcitEnergies.push_back(2715.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(97.0e-3*picosecond);

    ExcitEnergies.push_back(2771.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(107.0e-3*picosecond);

    ExcitEnergies.push_back(2908.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(3795.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(14.0*nanosecond);

    ExcitEnergies.push_back(4356.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(14.0*nanosecond);

}


G4Mg23GEMProbability::G4Mg23GEMProbability(const G4Mg23GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Mg23GEMProbability::copy_constructor meant to not be accessable");
}




const G4Mg23GEMProbability & G4Mg23GEMProbability::
operator=(const G4Mg23GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Mg23GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Mg23GEMProbability::operator==(const G4Mg23GEMProbability &) const
{
  return false;
}

G4bool G4Mg23GEMProbability::operator!=(const G4Mg23GEMProbability &) const
{
  return true;
}



