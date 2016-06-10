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
// $Id: G4Mg26GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Mg26GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Mg26GEMProbability::G4Mg26GEMProbability() :
  G4GEMProbability(26,12,0.0) // A,Z,Spin
{

    ExcitEnergies.push_back(1808.68*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(485.0e-3*picosecond);

    ExcitEnergies.push_back(2938.36*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(139.0e-3*picosecond);

    ExcitEnergies.push_back(3588.3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(6.6e-6*picosecond);

    ExcitEnergies.push_back(3940.5*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(6.2e-5*picosecond);

    ExcitEnergies.push_back(4318.4*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(201.0e-3*picosecond);

    ExcitEnergies.push_back(4331.9*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(55.0e-3*picosecond);

    ExcitEnergies.push_back(4349.8*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(104.0e-3*picosecond);

    ExcitEnergies.push_back(4834.3*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(35.0e-3*picosecond);

    ExcitEnergies.push_back(4900.3*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(49.0e-3*picosecond);

    ExcitEnergies.push_back(4972.2*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(444.0e-3*picosecond);

    ExcitEnergies.push_back(5290.8*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(69.0e-3*picosecond);

    ExcitEnergies.push_back(5473.9*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(24.0e-3*picosecond);

    ExcitEnergies.push_back(5690.1*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(49.0e-3*picosecond);

    ExcitEnergies.push_back(5715.5*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(94.0e-3*picosecond);

    ExcitEnergies.push_back(6256.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(52.0e-3*picosecond);

    ExcitEnergies.push_back(6621.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(69.0e-3*picosecond);

    ExcitEnergies.push_back(6744.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(55.0e-3*picosecond);

    ExcitEnergies.push_back(6877.7*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(83.0e-3*picosecond);

}

G4Mg26GEMProbability::~G4Mg26GEMProbability()
{}
