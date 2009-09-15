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
// $Id: G4Ne19GEMProbability.cc,v 1.6 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Ne19GEMProbability.hh"

G4Ne19GEMProbability::G4Ne19GEMProbability() :
  G4GEMProbability(19,10,1.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(238.33*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(18.0*nanosecond);

  ExcitEnergies.push_back(275.20*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(42.6*picosecond);

  ExcitEnergies.push_back(1507.8*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(2.8*picosecond);

  ExcitEnergies.push_back(1536.2*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(19.0e-3*picosecond);

  ExcitEnergies.push_back(1615.2*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(125.0e-3*picosecond);

  ExcitEnergies.push_back(2794.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(229.0e-3*picosecond);


  ExcitEnergies.push_back(10460.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(355.0*keV));

  ExcitEnergies.push_back(10480.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(45.0*keV));

  ExcitEnergies.push_back(11510.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(24.0*keV));

  ExcitEnergies.push_back(12230.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(200.0*keV));

  ExcitEnergies.push_back(12500.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(150.0*keV));

  ExcitEnergies.push_back(12690.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(180.0*keV));

}


G4Ne19GEMProbability::G4Ne19GEMProbability(const G4Ne19GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Ne19GEMProbability::copy_constructor meant to not be accessable");
}




const G4Ne19GEMProbability & G4Ne19GEMProbability::
operator=(const G4Ne19GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Ne19GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Ne19GEMProbability::operator==(const G4Ne19GEMProbability &) const
{
  return false;
}

G4bool G4Ne19GEMProbability::operator!=(const G4Ne19GEMProbability &) const
{
  return true;
}



