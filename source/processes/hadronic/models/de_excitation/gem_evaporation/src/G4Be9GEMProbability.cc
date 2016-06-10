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
// $Id: G4Be9GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Be9GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Be9GEMProbability::G4Be9GEMProbability() :
  G4GEMProbability(9,4,3.0/2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(1684.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(2429.4*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(0.77*keV));

  ExcitEnergies.push_back(2.78e3*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1080.0*keV));

  ExcitEnergies.push_back(3049.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(282.0*keV));

  ExcitEnergies.push_back(4704.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(743.0*keV));

  ExcitEnergies.push_back(6760.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1540.0*keV));

  ExcitEnergies.push_back(7940.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

  ExcitEnergies.push_back(11283.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(575.0*keV));

  ExcitEnergies.push_back(11810.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(400.0*keV));

  ExcitEnergies.push_back(13790.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(590.0*keV));

  ExcitEnergies.push_back(14392.2*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(0.381*keV));

  ExcitEnergies.push_back(14.4e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(800*keV));

  ExcitEnergies.push_back(15970.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(300.0*keV));

  ExcitEnergies.push_back(16671.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(41.0*keV));

  ExcitEnergies.push_back(16977.1*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(0.47*keV));

  ExcitEnergies.push_back(1730.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(17495.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(47.0*keV));

  ExcitEnergies.push_back(19200.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(310.0*keV));

  ExcitEnergies.push_back(20740.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));
}

G4Be9GEMProbability::~G4Be9GEMProbability()
{}
