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
// $Id: G4Ne23GEMProbability.cc,v 1.1 2002/06/06 18:02:33 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Ne23GEMProbability.hh"

G4Ne23GEMProbability::G4Ne23GEMProbability() :
  G4GEMProbability(23,10,5.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1017.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(178.0*picosecond);

  ExcitEnergies.push_back(1701.51*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(1822.5*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(2315.1*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(2517.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3221.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3431.8*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3458.2*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3830.9*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3836.8*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3988.2*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

}


G4Ne23GEMProbability::G4Ne23GEMProbability(const G4Ne23GEMProbability &right)
{
  G4Exception("G4Ne23GEMProbability::copy_constructor meant to not be accessable");
}




const G4Ne23GEMProbability & G4Ne23GEMProbability::
operator=(const G4Ne23GEMProbability &right)
{
  G4Exception("G4Ne23GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Ne23GEMProbability::operator==(const G4Ne23GEMProbability &right) const
{
  return false;
}

G4bool G4Ne23GEMProbability::operator!=(const G4Ne23GEMProbability &right) const
{
  return true;
}



