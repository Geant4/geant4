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
// $Id: G4Mg27GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Mg27GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Mg27GEMProbability::G4Mg27GEMProbability() :
  G4GEMProbability(27,12,1.0/2.0) // A,Z,Spin
{

    ExcitEnergies.push_back(984.66*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(1.4*picosecond);

    ExcitEnergies.push_back(1698.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(1.2*picosecond);

    ExcitEnergies.push_back(1940.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(1.1*picosecond);

    ExcitEnergies.push_back(3109.4*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(100.0e-3*picosecond);

    ExcitEnergies.push_back(3426.9*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(100.0e-3*picosecond);

    ExcitEnergies.push_back(3475.3*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(3490.7*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(15.0e-3*picosecond);

    ExcitEnergies.push_back(3559.2*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(3760.4*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(610.0e-3*picosecond);

    ExcitEnergies.push_back(3785.9*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(25.0e-3*picosecond);

    ExcitEnergies.push_back(3884.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(700.0e-3*picosecond);

    ExcitEnergies.push_back(4149.8*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(4398.2*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(65.0e-3*picosecond);

    ExcitEnergies.push_back(4552.8*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(15.0e-3*picosecond);

    ExcitEnergies.push_back(4827.3*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(4992.3*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(5028.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(40.0e-3*picosecond);

    ExcitEnergies.push_back(5172.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(15.0e-3*picosecond);

    ExcitEnergies.push_back(5372.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(25.0e-3*picosecond);

    ExcitEnergies.push_back(5422.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(5627.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(5764.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(25.0e-3*picosecond);

    ExcitEnergies.push_back(5821.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

}

G4Mg27GEMProbability::~G4Mg27GEMProbability() 
{}
