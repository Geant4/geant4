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
// $Id: G4O14GEMProbability.cc,v 1.4 2005/06/04 13:25:25 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4O14GEMProbability.hh"

G4O14GEMProbability::G4O14GEMProbability() :
  G4GEMProbability(14,8,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(5920.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(50.0*keV));

  ExcitEnergies.push_back(6272.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(103.0*keV));

  ExcitEnergies.push_back(6590.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(60.0*keV));

  ExcitEnergies.push_back(7768.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(76.0*keV));

  ExcitEnergies.push_back(9915.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

}


G4O14GEMProbability::G4O14GEMProbability(const G4O14GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O14GEMProbability::copy_constructor meant to not be accessable");
}




const G4O14GEMProbability & G4O14GEMProbability::
operator=(const G4O14GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O14GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4O14GEMProbability::operator==(const G4O14GEMProbability &) const
{
  return false;
}

G4bool G4O14GEMProbability::operator!=(const G4O14GEMProbability &) const
{
  return true;
}



