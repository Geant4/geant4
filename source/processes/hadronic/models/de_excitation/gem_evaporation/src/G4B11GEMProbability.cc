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
// $Id: G4B11GEMProbability.cc,v 1.6 2009-09-15 12:54:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4B11GEMProbability.hh"

G4B11GEMProbability::G4B11GEMProbability() :
  G4GEMProbability(11,5,3.0/2.0) // A,Z,Spin
{
    ExcitEnergies.push_back(2124.693*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(3.8e-15*s);

    ExcitEnergies.push_back(4444.89*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.82e-15*s);

    ExcitEnergies.push_back(5020.31*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(0.236e-15*s);

    ExcitEnergies.push_back(6742.9*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(15.0e-15*s);

    ExcitEnergies.push_back(6791.8*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(1.18e-15*s);

    ExcitEnergies.push_back(7285.51*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.4e-15*s);

    ExcitEnergies.push_back(7977.84*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(0.4e-15*s);

    ExcitEnergies.push_back(8560.3*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.49e-15*s);

    ExcitEnergies.push_back(8920.2*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(4.4e-15*s);

    ExcitEnergies.push_back(9185.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.9*eV));

    ExcitEnergies.push_back(9274.4*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(4.0*keV));

    ExcitEnergies.push_back(9876.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(110.0*keV));

    ExcitEnergies.push_back(10260.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(165.0*keV));

    ExcitEnergies.push_back(10330.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(110.0*keV));

    ExcitEnergies.push_back(10597.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

    ExcitEnergies.push_back(10960.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(4500.0*keV));

    ExcitEnergies.push_back(11265.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(110.0*keV));

    ExcitEnergies.push_back(11444.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(103.0*keV));

    ExcitEnergies.push_back(11886.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(200.0*keV));

    ExcitEnergies.push_back(12.0e3*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1000.0*keV));

    ExcitEnergies.push_back(12557.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(210.0*keV));

    ExcitEnergies.push_back(12916.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(155.0*keV));

    ExcitEnergies.push_back(13137.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(426.0*keV));

    ExcitEnergies.push_back(13.16e3*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(430.0*keV));

    ExcitEnergies.push_back(14.04e3*keV);
    ExcitSpins.push_back(11.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(500.0*keV));

    ExcitEnergies.push_back(14.34e3*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(254.0*keV));

    ExcitEnergies.push_back(14565.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30.0*keV));

    ExcitEnergies.push_back(15.32e3*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(635.0*keV));

    ExcitEnergies.push_back(16437.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30.0*keV));

    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4B11GEMProbability::G4B11GEMProbability(const G4B11GEMProbability &): G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4B11GEMProbability::copy_constructor meant to not be accessable");}




const G4B11GEMProbability & G4B11GEMProbability::
operator=(const G4B11GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4B11GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4B11GEMProbability::operator==(const G4B11GEMProbability &) const
{
  return false;
}

G4bool G4B11GEMProbability::operator!=(const G4B11GEMProbability &) const
{
  return true;
}



