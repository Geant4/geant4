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
// $Id: G4B13GEMProbability.cc,v 1.2 2003/05/30 13:23:23 hpw Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4B13GEMProbability.hh"

G4B13GEMProbability::G4B13GEMProbability() :
  G4GEMProbability(13,5,3.0/2.0) // A,Z,Spin
{

    ExcitEnergies.push_back(3534.7*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(0.2e-15*s);

    ExcitEnergies.push_back(3712*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(0.26e-15*s);

    ExcitEnergies.push_back(4131*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(0.04e-15*s);

    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4B13GEMProbability::G4B13GEMProbability(const G4B13GEMProbability &) : G4GEMProbability()
{
  G4Exception("G4B13GEMProbability::copy_constructor meant to not be accessable");}




const G4B13GEMProbability & G4B13GEMProbability::
operator=(const G4B13GEMProbability &)
{
  G4Exception("G4B13GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4B13GEMProbability::operator==(const G4B13GEMProbability &) const
{
  return false;
}

G4bool G4B13GEMProbability::operator!=(const G4B13GEMProbability &) const
{
  return true;
}



