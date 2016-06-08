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
// $Id: G4C15GEMProbability.cc,v 1.1 2002/06/06 17:59:56 larazb Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4C15GEMProbability.hh"

G4C15GEMProbability::G4C15GEMProbability() :
  G4GEMProbability(15,6,1.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(740.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(2.61e-9*s);

  ExcitEnergies.push_back(3105*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(40*keV));

  ExcitEnergies.push_back(4221*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(14*keV));

  ExcitEnergies.push_back(6370*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(20*keV));

  ExcitEnergies.push_back(6429*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(50*keV));

  ExcitEnergies.push_back(6461*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(14*keV));

  ExcitEnergies.push_back(6639*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(20*keV));

  ExcitEnergies.push_back(6845*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(14*keV));

  ExcitEnergies.push_back(6884*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(20*keV));

  ExcitEnergies.push_back(7098*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(15*keV));

  ExcitEnergies.push_back(7352*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(20*keV));

  ExcitEnergies.push_back(8470*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(40*keV));

  ExcitEnergies.push_back(8559*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(40*keV));

  ExcitEnergies.push_back(9789*keV);
  ExcitSpins.push_back(15.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(20*keV));

  ExcitEnergies.push_back(10248*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(20*keV));

  ExcitEnergies.push_back(11123*keV);
  ExcitSpins.push_back(19.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(30*keV));

  ExcitEnergies.push_back(11825*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(70*keV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4C15GEMProbability::G4C15GEMProbability(const G4C15GEMProbability &right)
{
  G4Exception("G4C15GEMProbability::copy_constructor meant to not be accessable");
}




const G4C15GEMProbability & G4C15GEMProbability::
operator=(const G4C15GEMProbability &right)
{
  G4Exception("G4C15GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4C15GEMProbability::operator==(const G4C15GEMProbability &right) const
{
  return false;
}

G4bool G4C15GEMProbability::operator!=(const G4C15GEMProbability &right) const
{
  return true;
}



