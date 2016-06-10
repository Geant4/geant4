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
// $Id: G4Ne22GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Ne22GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Ne22GEMProbability::G4Ne22GEMProbability() :
  G4GEMProbability(22,10,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1274.57*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(3.67*picosecond);

  ExcitEnergies.push_back(3357.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(225.0e-3*picosecond);

  ExcitEnergies.push_back(4456.7*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(17.0e-3*picosecond);

  ExcitEnergies.push_back(5147.5*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.8*picosecond);

  ExcitEnergies.push_back(5336.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(1.2e-3*picosecond);

  ExcitEnergies.push_back(5365.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(5523.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(35.0e-3*picosecond);

  ExcitEnergies.push_back(5641.3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(42.0e-3*picosecond);

  ExcitEnergies.push_back(5909.9*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(35.0e-3*picosecond);

  ExcitEnergies.push_back(6311.4*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(54.0e-3*picosecond);

  ExcitEnergies.push_back(6345.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(17.0e-3*picosecond);

  ExcitEnergies.push_back(6636.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(48.0e-3*picosecond);

  ExcitEnergies.push_back(6854.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(267.0e-6*picosecond);

  ExcitEnergies.push_back(7406.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(62.0e-3*picosecond);

  ExcitEnergies.push_back(423.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(47.0e-3*picosecond);

}

G4Ne22GEMProbability::~G4Ne22GEMProbability() 
{}
