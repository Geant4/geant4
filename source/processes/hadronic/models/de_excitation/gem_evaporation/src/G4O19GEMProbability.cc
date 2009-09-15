//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4O19GEMProbability.cc,v 1.6 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4O19GEMProbability.hh"

G4O19GEMProbability::G4O19GEMProbability() :
  G4GEMProbability(19,8,5.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(96.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(1.37e-3*picosecond);

  ExcitEnergies.push_back(1471.7*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(0.78*picosecond);

  ExcitEnergies.push_back(3154.5*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(1.0*picosecond);

  ExcitEnergies.push_back(4583.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(53.0*keV));

  ExcitEnergies.push_back(4707.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(15.0*keV));

  ExcitEnergies.push_back(5086.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(49.0*keV));

  ExcitEnergies.push_back(5149.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(3.4*keV));

  ExcitEnergies.push_back(5455.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(330.0*keV));

  ExcitEnergies.push_back(5706.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(7.8*keV));

  ExcitEnergies.push_back(6130.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(150.0*keV));

  ExcitEnergies.push_back(6200.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(140.0*keV));

  ExcitEnergies.push_back(6276.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(19.2*keV));
}


G4O19GEMProbability::G4O19GEMProbability(const G4O19GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O19GEMProbability::copy_constructor meant to not be accessable");
}




const G4O19GEMProbability & G4O19GEMProbability::
operator=(const G4O19GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O19GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4O19GEMProbability::operator==(const G4O19GEMProbability &) const
{
  return false;
}

G4bool G4O19GEMProbability::operator!=(const G4O19GEMProbability &) const
{
  return true;
}



