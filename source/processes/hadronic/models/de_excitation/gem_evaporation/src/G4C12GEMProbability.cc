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
// $Id: G4C12GEMProbability.cc,v 1.6 2009-09-15 12:54:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4C12GEMProbability.hh"

G4C12GEMProbability::G4C12GEMProbability() :
  G4GEMProbability(12,6,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(4438.91*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(10.8e-3*eV));

  ExcitEnergies.push_back(7654.2*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(8.5*keV));

  ExcitEnergies.push_back(9641*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(34*keV));

  ExcitEnergies.push_back(1.03E+4*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(3.0e3*keV));

  ExcitEnergies.push_back(10844*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(315*keV));

  ExcitEnergies.push_back(1.116E+4*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(430*keV));

  ExcitEnergies.push_back(1.1828E+4*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(260*keV));

  ExcitEnergies.push_back(12710*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(18.1*eV));

  ExcitEnergies.push_back(13352*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(375*keV));

  ExcitEnergies.push_back(14083*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(258*keV));

  ExcitEnergies.push_back(15110*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(42*eV));

  ExcitEnergies.push_back(16106.7*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(5.2*keV));

  ExcitEnergies.push_back(16.57E+3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300*keV));

  ExcitEnergies.push_back(17.23E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1150*keV));

  ExcitEnergies.push_back(17.76E+3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(80*keV));

  ExcitEnergies.push_back(18.13E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(600*keV));

  ExcitEnergies.push_back(18.27E+3*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300*keV));

  ExcitEnergies.push_back(18.38E+3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(220*keV));

  ExcitEnergies.push_back(18.39E+3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(43*keV));

  ExcitEnergies.push_back(18.6E+3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300*keV));

  ExcitEnergies.push_back(18.80E+3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100*keV));

  ExcitEnergies.push_back(19.2E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1100*keV));

  ExcitEnergies.push_back(19.39E+3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(800*keV));

  ExcitEnergies.push_back(19.55E+3*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(490*keV));

  ExcitEnergies.push_back(19.69E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(230*keV));

  ExcitEnergies.push_back(20.0E+3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100*keV));

  ExcitEnergies.push_back(20.27E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(140*keV));

  ExcitEnergies.push_back(20.5E+3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(180*keV));

  ExcitEnergies.push_back(20.62E+3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(200*keV));

  ExcitEnergies.push_back(21.60E+3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1200*keV));

  ExcitEnergies.push_back(22.0E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(800*keV));

  ExcitEnergies.push_back(22.40E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(275*keV));

  ExcitEnergies.push_back(22.65E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(3200*keV));

  ExcitEnergies.push_back(23.04E+3 *keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(60*keV));

  ExcitEnergies.push_back(23.52E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(230*keV));

  ExcitEnergies.push_back(23.92E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(400*keV));

  ExcitEnergies.push_back(25.30E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(510*keV));

  ExcitEnergies.push_back(25.4E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2000*keV));

  ExcitEnergies.push_back(27.0E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1400*keV));

  ExcitEnergies.push_back(27595.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30*keV));

  ExcitEnergies.push_back(28.2E+3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1600*keV));

  SetExcitationEnergiesPtr(&ExcitEnergies);
  SetExcitationSpinsPtr(&ExcitSpins);
  SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4C12GEMProbability::G4C12GEMProbability(const G4C12GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4C12GEMProbability::copy_constructor meant to not be accessable");
}




const G4C12GEMProbability & G4C12GEMProbability::
operator=(const G4C12GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4C12GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4C12GEMProbability::operator==(const G4C12GEMProbability &) const
{
  return false;
}

G4bool G4C12GEMProbability::operator!=(const G4C12GEMProbability &) const
{
  return true;
}



