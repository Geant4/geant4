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
// $Id: G4C13GEMProbability.cc,v 1.1 2002/06/06 17:59:54 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4C13GEMProbability.hh"

G4C13GEMProbability::G4C13GEMProbability() :
  G4GEMProbability(13,6,1.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(3088.4*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(1.04e-15*s);

  ExcitEnergies.push_back(3684.37*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(1.04e-15*s);

  ExcitEnergies.push_back( 3853.62*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(7.5e-12*s);

  ExcitEnergies.push_back(6864*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(6*keV));

  ExcitEnergies.push_back(7492*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(5*keV));

  ExcitEnergies.push_back(7547*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(1.2*keV));

  ExcitEnergies.push_back(7677*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(70*keV));

  ExcitEnergies.push_back(8.2E3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(1000*keV));

  ExcitEnergies.push_back(8860*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(150*keV));

  ExcitEnergies.push_back(9498*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(5*keV));

  ExcitEnergies.push_back(9897*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(26*keV));

  ExcitEnergies.push_back(10753*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(55*keV));

  ExcitEnergies.push_back(10818*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(24*keV));

  ExcitEnergies.push_back(10996*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(37*keV));

  ExcitEnergies.push_back(11080*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(4*keV));

  ExcitEnergies.push_back(11851*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(68*keV));

  ExcitEnergies.push_back(11970*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(200*keV));

  ExcitEnergies.push_back(12106*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(81*keV));

  ExcitEnergies.push_back(12400*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(150*keV));

  ExcitEnergies.push_back(13280*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(340*keV));

  ExcitEnergies.push_back(13410*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(35*keV));

  ExcitEnergies.push_back(13560*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(500*keV));

  ExcitEnergies.push_back(13760*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(300*keV));

  ExcitEnergies.push_back(14120*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(200*keV));

  ExcitEnergies.push_back(14.39E3*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(260*keV));

  ExcitEnergies.push_back(14940*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(380*keV));

  ExcitEnergies.push_back(15106*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(5*keV));

  ExcitEnergies.push_back(19500*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(450*keV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4C13GEMProbability::G4C13GEMProbability(const G4C13GEMProbability &right)
{
  G4Exception("G4C13GEMProbability::copy_constructor meant to not be accessable");
}




const G4C13GEMProbability & G4C13GEMProbability::
operator=(const G4C13GEMProbability &right)
{
  G4Exception("G4C13GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4C13GEMProbability::operator==(const G4C13GEMProbability &right) const
{
  return false;
}

G4bool G4C13GEMProbability::operator!=(const G4C13GEMProbability &right) const
{
  return true;
}



