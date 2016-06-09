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
// $Id: G4Be10GEMProbability.cc,v 1.1 2002/06/06 17:59:29 larazb Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Be10GEMProbability.hh"

G4Be10GEMProbability::G4Be10GEMProbability() :
  G4GEMProbability(10,4,0.0) // A,Z,Spin
{
  ExcitEnergies.push_back(3368.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(125.0e-15*s);

  ExcitEnergies.push_back(5958.3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(55.0e-15*s);

  ExcitEnergies.push_back(5959.9*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(0.8e-12*s);

  ExcitEnergies.push_back(7371.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(15.7*keV));

  ExcitEnergies.push_back(7542.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(6.3*keV));

  ExcitEnergies.push_back(9270.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(150.0*keV));

  ExcitEnergies.push_back(9400.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(291.0*keV));

  ExcitEnergies.push_back(11.76e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(121.0*keV));

  ExcitEnergies.push_back(17790.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(110.0*keV));

  ExcitEnergies.push_back(18550.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(350.0*keV));


  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Be10GEMProbability::G4Be10GEMProbability(const G4Be10GEMProbability &right)
{
  G4Exception("G4Be10GEMProbability::copy_constructor meant to not be accessable");
}




const G4Be10GEMProbability & G4Be10GEMProbability::
operator=(const G4Be10GEMProbability &right)
{
  G4Exception("G4Be10GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Be10GEMProbability::operator==(const G4Be10GEMProbability &right) const
{
  return false;
}

G4bool G4Be10GEMProbability::operator!=(const G4Be10GEMProbability &right) const
{
  return true;
}



