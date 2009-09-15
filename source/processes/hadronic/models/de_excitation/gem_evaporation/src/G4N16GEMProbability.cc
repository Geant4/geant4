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
// $Id: G4N16GEMProbability.cc,v 1.6 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4N16GEMProbability.hh"

G4N16GEMProbability::G4N16GEMProbability() :
  G4GEMProbability(16,7,2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(120.1*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(5.25e-6*s);

  ExcitEnergies.push_back(397.5*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(4.5e-12*s);

  ExcitEnergies.push_back( 3355*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(15*keV));

  ExcitEnergies.push_back( 3519*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(3*keV));

  ExcitEnergies.push_back( 3960*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2*keV));

  ExcitEnergies.push_back( 4319*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(20*keV));

  ExcitEnergies.push_back( 4387*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(82*keV));

  ExcitEnergies.push_back( 4760*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(250*keV));

  ExcitEnergies.push_back( 4776*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(59*keV));

  ExcitEnergies.push_back( 5050*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(19*keV));

  ExcitEnergies.push_back( 5130*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(7*keV));

  ExcitEnergies.push_back( 5150*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(7*keV));

  ExcitEnergies.push_back( 5232*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(4*keV));

  ExcitEnergies.push_back( 5240*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(260*keV));

  ExcitEnergies.push_back( 5250*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(320*keV));

  ExcitEnergies.push_back( 5518*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(7*keV));

  ExcitEnergies.push_back( 5730*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(7*keV));

  ExcitEnergies.push_back( 6009*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(270*keV));

  ExcitEnergies.push_back( 6168*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(7*keV));

  ExcitEnergies.push_back( 6373*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30*keV));

  ExcitEnergies.push_back( 6513*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(34*keV));

  ExcitEnergies.push_back( 6840*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(140*keV));

  ExcitEnergies.push_back( 7020*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(22*keV));

  ExcitEnergies.push_back( 7250*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(17*keV));

  ExcitEnergies.push_back( 7573*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(7*keV));

  ExcitEnergies.push_back( 7877*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100*keV));

  ExcitEnergies.push_back( 8365*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(18*keV));

  ExcitEnergies.push_back( 8490*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(50*keV));

  ExcitEnergies.push_back( 8720*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(40*keV));

  ExcitEnergies.push_back( 9160*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100*keV));

  ExcitEnergies.push_back( 9459*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100*keV));

  ExcitEnergies.push_back( 9928*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(12*keV));

  ExcitEnergies.push_back(10055*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30*keV));

  ExcitEnergies.push_back(10270*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(165*keV));

  ExcitEnergies.push_back(10710*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(120*keV));

  ExcitEnergies.push_back(11620*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(220*keV));

  ExcitEnergies.push_back(11701*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(12*keV));

  ExcitEnergies.push_back(14410*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(180*keV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4N16GEMProbability::G4N16GEMProbability(const G4N16GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4N16GEMProbability::copy_constructor meant to not be accessable");
}




const G4N16GEMProbability & G4N16GEMProbability::
operator=(const G4N16GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4N16GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4N16GEMProbability::operator==(const G4N16GEMProbability &) const
{
  return false;
}

G4bool G4N16GEMProbability::operator!=(const G4N16GEMProbability &) const
{
  return true;
}



