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
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Ne23GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Ne23GEMProbability::G4Ne23GEMProbability() :
  G4GEMProbability(23,10,5.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1017.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(178.0*picosecond);

  ExcitEnergies.push_back(1701.51*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(1822.5*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(2315.1*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(2517.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3221.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3431.8*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3458.2*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3830.9*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3836.8*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(3988.2*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

}

G4Ne23GEMProbability::~G4Ne23GEMProbability() 
{}

