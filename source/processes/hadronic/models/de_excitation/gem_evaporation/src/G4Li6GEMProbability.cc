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
// $Id: G4Li6GEMProbability.cc,v 1.6 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Li6GEMProbability.hh"

G4Li6GEMProbability::G4Li6GEMProbability() :
  G4GEMProbability(6,3,1.0) // A,Z,Spin
{
  ExcitEnergies.push_back(2186.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(24.0*keV));
  
  ExcitEnergies.push_back(3562.88*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(8.2*eV));

  ExcitEnergies.push_back(4310.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.7*MeV));

  ExcitEnergies.push_back(5366.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(540.0*keV));

  ExcitEnergies.push_back(5.65E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.5*MeV));

  ExcitEnergies.push_back(15800.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(17.8*MeV));

  ExcitEnergies.push_back(23.0E+3*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(12.0*MeV));

  ExcitEnergies.push_back(25.0E+3*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(4.0*MeV));
  
  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Li6GEMProbability::G4Li6GEMProbability(const G4Li6GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Li6GEMProbability::copy_constructor meant to not be accessable");
}




const G4Li6GEMProbability & G4Li6GEMProbability::
operator=(const G4Li6GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Li6GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Li6GEMProbability::operator==(const G4Li6GEMProbability &) const
{
  return false;
}

G4bool G4Li6GEMProbability::operator!=(const G4Li6GEMProbability &) const
{
  return true;
}



