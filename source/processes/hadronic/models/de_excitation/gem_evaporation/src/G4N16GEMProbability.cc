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
// $Id: G4N16GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4N16GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4N16GEMProbability::G4N16GEMProbability() :
  G4GEMProbability(16,7,2.0) // A,Z,Spin
{
  ExcitEnergies.push_back(120.42*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(5.25e-6*s);

  ExcitEnergies.push_back(298.22*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(9.13e-11*s);

  ExcitEnergies.push_back(397.27*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(4.5e-12*s);

  ExcitEnergies.push_back( 3355*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(15*keV));

  ExcitEnergies.push_back( 3519*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(3*keV));

  ExcitEnergies.push_back( 3960*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(2*keV));

  ExcitEnergies.push_back( 4319*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back( 4387*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(82*keV));

  ExcitEnergies.push_back( 4760*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(250*keV));

  ExcitEnergies.push_back( 4776*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(59*keV));

  ExcitEnergies.push_back( 5050*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(19*keV));

  ExcitEnergies.push_back( 5130*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(7*keV));

  ExcitEnergies.push_back( 5150*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(7*keV));

  ExcitEnergies.push_back( 5232*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(4*keV));

  ExcitEnergies.push_back( 5240*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(260*keV));

  ExcitEnergies.push_back( 5250*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(320*keV));

  ExcitEnergies.push_back( 5518*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(7*keV));

  ExcitEnergies.push_back( 5730*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(7*keV));

  ExcitEnergies.push_back( 6009*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(270*keV));

  ExcitEnergies.push_back( 6168*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(7*keV));

  ExcitEnergies.push_back( 6373*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));

  ExcitEnergies.push_back( 6513*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(34*keV));

  ExcitEnergies.push_back( 6840*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(140*keV));

  ExcitEnergies.push_back( 7020*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(22*keV));

  ExcitEnergies.push_back( 7250*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(17*keV));

  ExcitEnergies.push_back( 7573*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(7*keV));

  ExcitEnergies.push_back( 7877*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(100*keV));

  ExcitEnergies.push_back( 8365*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(18*keV));

  ExcitEnergies.push_back( 8490*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(50*keV));

  ExcitEnergies.push_back( 8720*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(40*keV));

  ExcitEnergies.push_back( 9160*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(100*keV));

  ExcitEnergies.push_back( 9459*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(100*keV));

  ExcitEnergies.push_back( 9928*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(12*keV));

  ExcitEnergies.push_back(10055*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));

  ExcitEnergies.push_back(10270*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(165*keV));

  ExcitEnergies.push_back(10710*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(120*keV));

  ExcitEnergies.push_back(11620*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(220*keV));

  ExcitEnergies.push_back(11701*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(12*keV));

  ExcitEnergies.push_back(14410*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(180*keV));
}

G4N16GEMProbability::~G4N16GEMProbability() 
{}
