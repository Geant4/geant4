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
// $Id: G4Be10GEMProbability.cc,v 1.6 2009-09-15 12:54:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(15.7*keV));

  ExcitEnergies.push_back(7542.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(6.3*keV));

  ExcitEnergies.push_back(9270.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(150.0*keV));

  ExcitEnergies.push_back(9400.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(291.0*keV));

  ExcitEnergies.push_back(11.76e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(121.0*keV));

  ExcitEnergies.push_back(17790.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(110.0*keV));

  ExcitEnergies.push_back(18550.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(350.0*keV));


  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4Be10GEMProbability::G4Be10GEMProbability(const G4Be10GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Be10GEMProbability::copy_constructor meant to not be accessable");
}




const G4Be10GEMProbability & G4Be10GEMProbability::
operator=(const G4Be10GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Be10GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Be10GEMProbability::operator==(const G4Be10GEMProbability &) const
{
  return false;
}

G4bool G4Be10GEMProbability::operator!=(const G4Be10GEMProbability &) const
{
  return true;
}



