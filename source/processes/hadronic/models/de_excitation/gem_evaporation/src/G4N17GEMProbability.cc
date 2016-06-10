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
// $Id: G4N17GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4N17GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4N17GEMProbability::G4N17GEMProbability() :
  G4GEMProbability(17,7,1.0/2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(1373.9*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(64e-15*s);

  ExcitEnergies.push_back(1849.6*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(28e-12*s);

  ExcitEnergies.push_back(1906.8*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(7.60e-12*s);

  ExcitEnergies.push_back(2526.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(23e-12*s);

  ExcitEnergies.push_back(3128.9*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(191e-15*s);

  ExcitEnergies.push_back(3204.2*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(21e-15*s);

  ExcitEnergies.push_back(3628.7*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(8.30e-12*s);

  ExcitEnergies.push_back(3663.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(243e-15*s);

  ExcitEnergies.push_back(3906.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(36e-15*s);

  ExcitEnergies.push_back(4006.4*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(10e-15*s);

  ExcitEnergies.push_back(4208.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(49e-15*s);

  ExcitEnergies.push_back(4415.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(42e-15*s);

  ExcitEnergies.push_back(5170.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(42e-15*s);

  ExcitEnergies.push_back(5195.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(66e-15*s);

  ExcitEnergies.push_back(5514.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69e-15*s);

  ExcitEnergies.push_back(5770.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(83e-15*s);
}

G4N17GEMProbability::~G4N17GEMProbability() 
{}
