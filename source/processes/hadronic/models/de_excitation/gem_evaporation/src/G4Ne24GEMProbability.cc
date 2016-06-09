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
// $Id: G4Ne24GEMProbability.cc,v 1.2 2003/11/03 17:53:04 hpw Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Ne24GEMProbability.hh"

G4Ne24GEMProbability::G4Ne24GEMProbability() :
  G4GEMProbability(24,10,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1980.8*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(658.0e-3*picosecond);

  ExcitEnergies.push_back(3867.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3962.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(21.0*nanosecond);

  ExcitEnergies.push_back(4764.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(2.3*picosecond);

  ExcitEnergies.push_back(5576.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(21.0*nanosecond);
}


G4Ne24GEMProbability::G4Ne24GEMProbability(const G4Ne24GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Ne24GEMProbability::copy_constructor meant to not be accessable");
}




const G4Ne24GEMProbability & G4Ne24GEMProbability::
operator=(const G4Ne24GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Ne24GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Ne24GEMProbability::operator==(const G4Ne24GEMProbability &) const
{
  return false;
}

G4bool G4Ne24GEMProbability::operator!=(const G4Ne24GEMProbability &) const
{
  return true;
}



