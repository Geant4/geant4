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
// $Id: G4O15GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4O15GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4O15GEMProbability::G4O15GEMProbability() :
  G4GEMProbability(15,8,1.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(5183.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(0.07*picosecond);

  ExcitEnergies.push_back(5240.9*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(2.2*picosecond);

  ExcitEnergies.push_back(6176.3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(33.0E-3*picosecond);

  ExcitEnergies.push_back(6793.1*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(19.0E-3*picosecond);

  ExcitEnergies.push_back(6859.4*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(0.07*picosecond);

  ExcitEnergies.push_back(7556.8*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.6*keV));

  ExcitEnergies.push_back(8284.3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.6*keV));

  ExcitEnergies.push_back(8743.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(32.0*keV));

  ExcitEnergies.push_back(8922.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(8927.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(8982.4*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.9*keV));

  ExcitEnergies.push_back(9487.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(10.1*keV));

  ExcitEnergies.push_back(9527.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(280.0*keV));

  ExcitEnergies.push_back(9610.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(8.8*keV));

  ExcitEnergies.push_back(9662.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2.0*keV));

  ExcitEnergies.push_back(9720.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1185.0*keV));

  ExcitEnergies.push_back(10291.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.0*keV));

  ExcitEnergies.push_back(10296.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(11.0*keV));

  ExcitEnergies.push_back(10478.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(25.0*keV));

  ExcitEnergies.push_back(10506.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(140.0*keV));

  ExcitEnergies.push_back(10917.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(90.0*keV));

  ExcitEnergies.push_back(10938.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(99.0*keV));

  ExcitEnergies.push_back(11025.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(25.0*keV));

  ExcitEnergies.push_back(11218.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(11569.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20.0*keV));

  ExcitEnergies.push_back(11616.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(11748.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(99.0*keV));

  ExcitEnergies.push_back(11846.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(65.0*keV));

  ExcitEnergies.push_back(11980.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20.0*keV));

  ExcitEnergies.push_back(12129.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(12471.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(77.0*keV));

  ExcitEnergies.push_back(12835.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(16.0*keV));

  ExcitEnergies.push_back(13450.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

  ExcitEnergies.push_back(14030.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(160.0*keV));

  ExcitEnergies.push_back(14270.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(340.0*keV));

  ExcitEnergies.push_back(14465.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(15100.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

  ExcitEnergies.push_back(15900.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(350.0*keV));

  ExcitEnergies.push_back(16430.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(560.0*keV));

  ExcitEnergies.push_back(16900.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

  ExcitEnergies.push_back(17510.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(600.0*keV));

  ExcitEnergies.push_back(17990.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(18400.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

  ExcitEnergies.push_back(20500.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2000.0*keV));

  ExcitEnergies.push_back(22000.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2000.0*keV));

  ExcitEnergies.push_back(26000.0*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(600.0*keV));

  ExcitEnergies.push_back(28000.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2500.0*keV));
}

G4O15GEMProbability::~G4O15GEMProbability()
{}
