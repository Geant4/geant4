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
// $Id: G4B10GEMProbability.cc,v 1.6 2009-09-15 12:54:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4B10GEMProbability.hh"

G4B10GEMProbability::G4B10GEMProbability() :
  G4GEMProbability(10,5,3.0) // A,Z,Spin
{
    ExcitEnergies.push_back(718.35*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(0.707e-9*s);

    ExcitEnergies.push_back(1740.15*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(5.0e-15*s);

    ExcitEnergies.push_back(2154.3*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(1.48e-12*s);

    ExcitEnergies.push_back(3587.1*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(106.0e-15*s);

    ExcitEnergies.push_back(4774.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(8.7*keV));

    ExcitEnergies.push_back(5110.3*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.98*keV));

    ExcitEnergies.push_back(5163.9*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(4.2e-15*s);

    ExcitEnergies.push_back(5180.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(110.0*keV));

    ExcitEnergies.push_back(5919.5*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(6.0*keV));

    ExcitEnergies.push_back(6025.0*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.05*keV));

    ExcitEnergies.push_back(6127.2*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.36*keV));

    ExcitEnergies.push_back(6561.0*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(25.1*keV));

    ExcitEnergies.push_back(6873.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(120.0*keV));

    ExcitEnergies.push_back(7002.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

    ExcitEnergies.push_back(7430.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

    ExcitEnergies.push_back(7467.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(65.0*keV));

    ExcitEnergies.push_back(7479.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(74.0*keV));

    ExcitEnergies.push_back(7560.7*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.65*keV));

    ExcitEnergies.push_back(7670.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(250.0*keV));

    ExcitEnergies.push_back(7819.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(260.0*keV));

    ExcitEnergies.push_back(8070.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(800.0*keV));

    ExcitEnergies.push_back(8700.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(200.0*keV));

    ExcitEnergies.push_back(8889.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(84.0*keV));

    ExcitEnergies.push_back(8895.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(40.0*keV));

    ExcitEnergies.push_back(9700.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(700.0*keV));

    ExcitEnergies.push_back(10840.0*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

    ExcitEnergies.push_back(11520.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(500.0*keV));

    ExcitEnergies.push_back(12560.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

    ExcitEnergies.push_back(13490.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

    ExcitEnergies.push_back(14.4e3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(800.0*keV));

    ExcitEnergies.push_back(18.2e3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1500.0*keV));

    ExcitEnergies.push_back(18430.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(340.0*keV));

    ExcitEnergies.push_back(18800.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(600.0*keV));

    ExcitEnergies.push_back(19290.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(190.0*keV));


    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4B10GEMProbability::G4B10GEMProbability(const G4B10GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4B10GEMProbability::copy_constructor meant to not be accessable");}




const G4B10GEMProbability & G4B10GEMProbability::
operator=(const G4B10GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4B10GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4B10GEMProbability::operator==(const G4B10GEMProbability &) const
{
  return false;
}

G4bool G4B10GEMProbability::operator!=(const G4B10GEMProbability &) const
{
  return true;
}



