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
// $Id: G4Li7GEMProbability.cc,v 1.3 2004/12/07 13:47:13 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Li7GEMProbability.hh"

G4Li7GEMProbability::G4Li7GEMProbability() :
  G4GEMProbability(7,3,3.0/2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(477.612*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(73.0E-15*second);
  
  ExcitEnergies.push_back(4630.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(93.0*keV));

  ExcitEnergies.push_back(6680.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.9*MeV));

  ExcitEnergies.push_back(7459.7*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(89.0*keV));

  ExcitEnergies.push_back(9.67E+3*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(400.0*keV));

  ExcitEnergies.push_back(98500.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1200.0*keV));

  ExcitEnergies.push_back(11240.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(260.0*keV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Li7GEMProbability::G4Li7GEMProbability(const G4Li7GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Li7GEMProbability::copy_constructor meant to not be accessable");
}




const G4Li7GEMProbability & G4Li7GEMProbability::
operator=(const G4Li7GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Li7GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Li7GEMProbability::operator==(const G4Li7GEMProbability &) const
{
  return false;
}

G4bool G4Li7GEMProbability::operator!=(const G4Li7GEMProbability &) const
{
  return true;
}



