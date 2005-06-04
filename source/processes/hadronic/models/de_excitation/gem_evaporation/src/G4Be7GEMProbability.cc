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
// $Id: G4Be7GEMProbability.cc,v 1.4 2005-06-04 13:25:25 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Be7GEMProbability.hh"

G4Be7GEMProbability::G4Be7GEMProbability() :
  G4GEMProbability(7,4,3.0/2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(429.08*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(133.0e-15*s);

  ExcitEnergies.push_back(4570.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(175.0*keV));

  ExcitEnergies.push_back(6.73e3*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.2*MeV));

  ExcitEnergies.push_back(7210.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.5*MeV));

  ExcitEnergies.push_back(9900.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.8*MeV));

  ExcitEnergies.push_back(11010.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(320.0*keV));

  ExcitEnergies.push_back(17000.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(6.5*MeV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Be7GEMProbability::G4Be7GEMProbability(const G4Be7GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Be7GEMProbability::copy_constructor meant to not be accessable");
}




const G4Be7GEMProbability & G4Be7GEMProbability::
operator=(const G4Be7GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Be7GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Be7GEMProbability::operator==(const G4Be7GEMProbability &) const
{
  return false;
}

G4bool G4Be7GEMProbability::operator!=(const G4Be7GEMProbability &) const
{
  return true;
}



