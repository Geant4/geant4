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
// $Id: G4N12GEMProbability.cc,v 1.2 2003/11/03 17:53:04 hpw Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4N12GEMProbability.hh"

G4N12GEMProbability::G4N12GEMProbability() :
  G4GEMProbability(12,7,1.0) // A,Z,Spin
{

  ExcitEnergies.push_back(960.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(20*keV));

  ExcitEnergies.push_back(1189*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(100*keV));

  ExcitEnergies.push_back(2415*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(45*keV));

  ExcitEnergies.push_back(3118*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(210*keV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4N12GEMProbability::G4N12GEMProbability(const G4N12GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4N12GEMProbability::copy_constructor meant to not be accessable");
}




const G4N12GEMProbability & G4N12GEMProbability::
operator=(const G4N12GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4N12GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4N12GEMProbability::operator==(const G4N12GEMProbability &) const
{
  return false;
}

G4bool G4N12GEMProbability::operator!=(const G4N12GEMProbability &) const
{
  return true;
}



