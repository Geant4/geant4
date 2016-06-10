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
// $Id: G4Mg28GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Mg28GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Mg28GEMProbability::G4Mg28GEMProbability() :
  G4GEMProbability(28,12,0.0) // A,Z,Spin
{

    ExcitEnergies.push_back(1473.4*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(1.7*picosecond);

    ExcitEnergies.push_back(3862.7*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(800.0e-3*picosecond);

    ExcitEnergies.push_back(4020.2*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(150.0e-3*picosecond);

    ExcitEnergies.push_back(4557.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(40.0e-3*picosecond);

    ExcitEnergies.push_back(4878.6*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(120.0e-3*picosecond);

    ExcitEnergies.push_back(5171.8*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(170.0e-3*picosecond);

    ExcitEnergies.push_back(5192.7*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(30.0e-3*picosecond);

    ExcitEnergies.push_back(5271.7*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(150.0e-3*picosecond);

    ExcitEnergies.push_back(5702.3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(300.0e-3*picosecond);

}

G4Mg28GEMProbability::~G4Mg28GEMProbability() 
{}
