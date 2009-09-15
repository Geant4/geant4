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
// $Id: G4Mg23GEMProbability.cc,v 1.5 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Mg23GEMProbability.hh"

G4Mg23GEMProbability::G4Mg23GEMProbability() :
  G4GEMProbability(23,12,3.0/2.0) // A,Z,Spin
{

    ExcitEnergies.push_back(450.70*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(1.25*picosecond);

    ExcitEnergies.push_back(2051.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(55.0e-3*picosecond);

    ExcitEnergies.push_back(2359.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(575.0e-3*picosecond);

    ExcitEnergies.push_back(2715.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(97.0e-3*picosecond);

    ExcitEnergies.push_back(2771.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(107.0e-3*picosecond);

    ExcitEnergies.push_back(2908.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(3795.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(14.0*nanosecond);

    ExcitEnergies.push_back(4356.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(14.0*nanosecond);

}


G4Mg23GEMProbability::G4Mg23GEMProbability(const G4Mg23GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Mg23GEMProbability::copy_constructor meant to not be accessable");
}




const G4Mg23GEMProbability & G4Mg23GEMProbability::
operator=(const G4Mg23GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Mg23GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Mg23GEMProbability::operator==(const G4Mg23GEMProbability &) const
{
  return false;
}

G4bool G4Mg23GEMProbability::operator!=(const G4Mg23GEMProbability &) const
{
  return true;
}



