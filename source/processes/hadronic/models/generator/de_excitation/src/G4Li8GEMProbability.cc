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
// $Id: G4Li8GEMProbability.cc,v 1.1 2002/06/06 18:02:07 larazb Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Li8GEMProbability.hh"

G4Li8GEMProbability::G4Li8GEMProbability() :
  G4GEMProbability(8,3,2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(980.8*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(8.0E-15*second);
  
  ExcitEnergies.push_back(2255.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(33.0*keV));

  ExcitEnergies.push_back(3210.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(1.0*MeV));

  ExcitEnergies.push_back(5400.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(650.0*keV));

  ExcitEnergies.push_back(6.1e3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(1.0*MeV));

  ExcitEnergies.push_back(6530.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(35.0*keV));

  ExcitEnergies.push_back(7.1e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(400.0*keV));

  ExcitEnergies.push_back(9000.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(6000.0*keV));

  ExcitEnergies.push_back(10822.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(12.0*keV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Li8GEMProbability::G4Li8GEMProbability(const G4Li8GEMProbability &right)
{
  G4Exception("G4Li8GEMProbability::copy_constructor meant to not be accessable");
}




const G4Li8GEMProbability & G4Li8GEMProbability::
operator=(const G4Li8GEMProbability &right)
{
  G4Exception("G4Li8GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Li8GEMProbability::operator==(const G4Li8GEMProbability &right) const
{
  return false;
}

G4bool G4Li8GEMProbability::operator!=(const G4Li8GEMProbability &right) const
{
  return true;
}



