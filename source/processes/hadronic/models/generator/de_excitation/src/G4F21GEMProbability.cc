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
// $Id: G4F21GEMProbability.cc,v 1.1 2002/06/06 18:00:53 larazb Exp $
// GEANT4 tag $Name: geant4-05-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4F21GEMProbability.hh"

G4F21GEMProbability::G4F21GEMProbability() :
  G4GEMProbability(21,9,5.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(279.92*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(6.1*nanosecond);

  ExcitEnergies.push_back(1100.9*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(305.0e-3*picosecond);

  ExcitEnergies.push_back(3449.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(0.7*picosecond);

  ExcitEnergies.push_back(3508.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(0.7*picosecond);

}


G4F21GEMProbability::G4F21GEMProbability(const G4F21GEMProbability &right)
{
  G4Exception("G4F21GEMProbability::copy_constructor meant to not be accessable");
}




const G4F21GEMProbability & G4F21GEMProbability::
operator=(const G4F21GEMProbability &right)
{
  G4Exception("G4F21GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4F21GEMProbability::operator==(const G4F21GEMProbability &right) const
{
  return false;
}

G4bool G4F21GEMProbability::operator!=(const G4F21GEMProbability &right) const
{
  return true;
}



