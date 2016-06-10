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
// $Id: G4O19GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4O19GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4O19GEMProbability::G4O19GEMProbability() :
  G4GEMProbability(19,8,5.0/2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(96.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(1.37e-3*picosecond);

  ExcitEnergies.push_back(1471.7*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(0.78*picosecond);

  ExcitEnergies.push_back(3154.5*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(1.0*picosecond);

  ExcitEnergies.push_back(4583.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(53.0*keV));

  ExcitEnergies.push_back(4707.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(15.0*keV));

  ExcitEnergies.push_back(5086.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(49.0*keV));

  ExcitEnergies.push_back(5149.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.4*keV));

  ExcitEnergies.push_back(5455.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(330.0*keV));

  ExcitEnergies.push_back(5706.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(7.8*keV));

  ExcitEnergies.push_back(6130.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(6200.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(140.0*keV));

  ExcitEnergies.push_back(6276.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(19.2*keV));
}

G4O19GEMProbability::~G4O19GEMProbability() 
{}

