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
// $Id: G4Be12GEMProbability.cc,v 1.2 2003/05/30 13:23:23 hpw Exp $
// GEANT4 tag $Name: geant4-05-02 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Be12GEMProbability.hh"

G4Be12GEMProbability::G4Be12GEMProbability() :
  G4GEMProbability(9,4,0.0) // A,Z,Spin
{
  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Be12GEMProbability::G4Be12GEMProbability(const G4Be12GEMProbability &) : G4GEMProbability()
{
  G4Exception("G4Be12GEMProbability::copy_constructor meant to not be accessable");
}




const G4Be12GEMProbability & G4Be12GEMProbability::
operator=(const G4Be12GEMProbability &)
{
  G4Exception("G4Be12GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Be12GEMProbability::operator==(const G4Be12GEMProbability &) const
{
  return false;
}

G4bool G4Be12GEMProbability::operator!=(const G4Be12GEMProbability &) const
{
  return true;
}



