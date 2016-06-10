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
// $Id: G4F20GEMProbability.cc 74869 2013-10-23 09:26:17Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4F20GEMProbability.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"

G4F20GEMProbability::G4F20GEMProbability() :
  G4GEMProbability(20,9,2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(655.95*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(0.28*picosecond);

  ExcitEnergies.push_back(822.9*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(55*picosecond);

  ExcitEnergies.push_back(983.8*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(1.1*picosecond);

  ExcitEnergies.push_back(1056.93*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(31.0e-3*picosecond);

  ExcitEnergies.push_back(1309.22*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.62*picosecond);

  ExcitEnergies.push_back(1843.4*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(2043.9*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(26.0e-3*picosecond);

  ExcitEnergies.push_back(2194.6*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(8.0e-3*picosecond);

  ExcitEnergies.push_back(2966.2*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(42.0e-3*picosecond);

  ExcitEnergies.push_back(3488.4*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(30.0e-3*picosecond);

  ExcitEnergies.push_back(3525.9*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(3587.1*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(6627.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*G4Log(2.0)/(0.29*keV));

  ExcitEnergies.push_back(6648.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(1.62*keV));

  ExcitEnergies.push_back(6685.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(3.80*keV));

  ExcitEnergies.push_back(6692.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(5.23*keV));

  ExcitEnergies.push_back(6696.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(1.05*keV));

  ExcitEnergies.push_back(6699.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(2.85*keV));

  ExcitEnergies.push_back(6709.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(1.14*keV));

  ExcitEnergies.push_back(6717.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(0.95*keV));

  ExcitEnergies.push_back(6791.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(1.9*keV));

  ExcitEnergies.push_back(6835.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(1.7*keV));

  ExcitEnergies.push_back(6837.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(0.4*keV));

  ExcitEnergies.push_back(6856.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(1.3*keV));

  ExcitEnergies.push_back(6858.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(19.0*keV));

  ExcitEnergies.push_back(7005.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(24.0*keV));

  ExcitEnergies.push_back(7076.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(24.0*keV));

  ExcitEnergies.push_back(7171.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(14.0*keV));

  ExcitEnergies.push_back(7311.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(33.0*keV));

  ExcitEnergies.push_back(7355.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(19.0*keV));

  ExcitEnergies.push_back(7410.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(10.0*keV));

  ExcitEnergies.push_back(7489.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(57.0*keV));

  ExcitEnergies.push_back(7503.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(85.0*keV));

  ExcitEnergies.push_back(7670.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(7800.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(8150.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(190.0*keV));

  ExcitEnergies.push_back(10228.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(10641.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(10807.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(330.0*keV));

}

G4F20GEMProbability::~G4F20GEMProbability()
{}

