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
// $Id: G4Li8GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Li8GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Li8GEMProbability::G4Li8GEMProbability() :
  G4GEMProbability(8,3,2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(980.8*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(8.0E-15*second);
  
  ExcitEnergies.push_back(2255.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(33.0*keV));

  ExcitEnergies.push_back(3210.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(1.0*MeV));

  ExcitEnergies.push_back(5400.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(650.0*keV));

  ExcitEnergies.push_back(6.1e3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(1.0*MeV));

  ExcitEnergies.push_back(6530.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(35.0*keV));

  ExcitEnergies.push_back(7.1e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(400.0*keV));

  ExcitEnergies.push_back(9000.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(6000.0*keV));

  ExcitEnergies.push_back(10822.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(12.0*keV));
}

G4Li8GEMProbability::~G4Li8GEMProbability()
{}
