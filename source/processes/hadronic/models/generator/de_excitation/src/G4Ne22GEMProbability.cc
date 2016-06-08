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
// $Id: G4Ne22GEMProbability.cc,v 1.1 2002/06/06 18:02:32 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Ne22GEMProbability.hh"

G4Ne22GEMProbability::G4Ne22GEMProbability() :
  G4GEMProbability(22,10,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1274.57*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(3.67*picosecond);

  ExcitEnergies.push_back(3357.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(225.0e-3*picosecond);

  ExcitEnergies.push_back(4456.7*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(17.0e-3*picosecond);

  ExcitEnergies.push_back(5147.5*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.8*picosecond);

  ExcitEnergies.push_back(5336.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(1.2e-3*picosecond);

  ExcitEnergies.push_back(5365.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(5523.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(35.0e-3*picosecond);

  ExcitEnergies.push_back(5641.3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(42.0e-3*picosecond);

  ExcitEnergies.push_back(5909.9*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(35.0e-3*picosecond);

  ExcitEnergies.push_back(6311.4*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(54.0e-3*picosecond);

  ExcitEnergies.push_back(6345.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(17.0e-3*picosecond);

  ExcitEnergies.push_back(6636.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(48.0e-3*picosecond);

  ExcitEnergies.push_back(6854.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(267.0e-6*picosecond);

  ExcitEnergies.push_back(7406.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(62.0e-3*picosecond);

  ExcitEnergies.push_back(423.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(47.0e-3*picosecond);

}


G4Ne22GEMProbability::G4Ne22GEMProbability(const G4Ne22GEMProbability &right)
{
  G4Exception("G4Ne22GEMProbability::copy_constructor meant to not be accessable");
}




const G4Ne22GEMProbability & G4Ne22GEMProbability::
operator=(const G4Ne22GEMProbability &right)
{
  G4Exception("G4Ne22GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Ne22GEMProbability::operator==(const G4Ne22GEMProbability &right) const
{
  return false;
}

G4bool G4Ne22GEMProbability::operator!=(const G4Ne22GEMProbability &right) const
{
  return true;
}



