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
// $Id: G4Be11GEMProbability.cc,v 1.3 2004/12/07 13:46:55 gunter Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Be11GEMProbability.hh"

G4Be11GEMProbability::G4Be11GEMProbability() :
  G4GEMProbability(11,4,1.0/2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(320.04*keV);
  ExcitSpins.push_back(1.0/2);
  ExcitLifetimes.push_back(115.0e-15*s);

  ExcitEnergies.push_back(1778.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

  ExcitEnergies.push_back(2.69e3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(200.0*keV));

  ExcitEnergies.push_back(3.41e3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(125.0*keV));

  ExcitEnergies.push_back(3887.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(10.0*keV));

  ExcitEnergies.push_back(3956.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(15.0*keV));

  ExcitEnergies.push_back(5240.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(45.0*keV));

  ExcitEnergies.push_back(5.86e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

  ExcitEnergies.push_back(6.51e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(120.0*keV));

  ExcitEnergies.push_back(6705.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(40.0*keV));

  ExcitEnergies.push_back(7.03e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

  ExcitEnergies.push_back(8816.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

  ExcitEnergies.push_back(10.59e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(210.0*keV));


  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Be11GEMProbability::G4Be11GEMProbability(const G4Be11GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Be11GEMProbability::copy_constructor meant to not be accessable");
}




const G4Be11GEMProbability & G4Be11GEMProbability::
operator=(const G4Be11GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Be11GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Be11GEMProbability::operator==(const G4Be11GEMProbability &) const
{
  return false;
}

G4bool G4Be11GEMProbability::operator!=(const G4Be11GEMProbability &) const
{
  return true;
}



