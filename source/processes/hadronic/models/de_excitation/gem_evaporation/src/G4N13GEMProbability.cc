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
// $Id: G4N13GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4N13GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4N13GEMProbability::G4N13GEMProbability() :
  G4GEMProbability(13,7,1.0/2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(2364.9*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(33*keV));

  ExcitEnergies.push_back(3502*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(66*keV));

  ExcitEnergies.push_back(3547*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(47*keV));

  ExcitEnergies.push_back(6364*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(11*keV));

  ExcitEnergies.push_back(6885*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(115*keV));

  ExcitEnergies.push_back(7155*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(9.0*keV));

  ExcitEnergies.push_back(7376*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(75*keV));

  ExcitEnergies.push_back(8000*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1500*keV));

  ExcitEnergies.push_back(8918*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(230*keV));

  ExcitEnergies.push_back(9476*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));

  ExcitEnergies.push_back(10.25E3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(270*keV));

  ExcitEnergies.push_back(10360*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));

  ExcitEnergies.push_back(10360*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(76*keV));

  ExcitEnergies.push_back(11530*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(430*keV));

  ExcitEnergies.push_back(11740*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(115*keV));

  ExcitEnergies.push_back(11740*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(250*keV));

  ExcitEnergies.push_back(11740*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(530*keV));

  ExcitEnergies.push_back(11878*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(380*keV));

  ExcitEnergies.push_back(12130*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(250*keV));

  ExcitEnergies.push_back(13.5E3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6500*keV));

  ExcitEnergies.push_back(14050*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(155*keV));

  ExcitEnergies.push_back(15064.6*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.10*keV));

  ExcitEnergies.push_back(15.3E3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(350*keV));

  ExcitEnergies.push_back(16020*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(130*keV));

  ExcitEnergies.push_back(18230*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(300*keV));

  ExcitEnergies.push_back(18406*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(66*keV));

  ExcitEnergies.push_back(18961*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(23*keV));

  ExcitEnergies.push_back(19830*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000*keV));

  ExcitEnergies.push_back(19880*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(520*keV));

  ExcitEnergies.push_back(20.9E3*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1100*keV));

  ExcitEnergies.push_back(22.4E3*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2000*keV));
}

G4N13GEMProbability::~G4N13GEMProbability()
{}
