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
// $Id: G4C15GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4C15GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4C15GEMProbability::G4C15GEMProbability() :
  G4GEMProbability(15,6,1.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(740.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(2.61e-9*s);

  ExcitEnergies.push_back(3105*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(40*keV));

  ExcitEnergies.push_back(4221*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(14*keV));

  ExcitEnergies.push_back(6370*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back(6429*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(50*keV));

  ExcitEnergies.push_back(6461*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(14*keV));

  ExcitEnergies.push_back(6639*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back(6845*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(14*keV));

  ExcitEnergies.push_back(6884*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back(7098*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(15*keV));

  ExcitEnergies.push_back(7352*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back(8470*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(40*keV));

  ExcitEnergies.push_back(8559*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(40*keV));

  ExcitEnergies.push_back(9789*keV);
  ExcitSpins.push_back(15.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back(10248*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back(11123*keV);
  ExcitSpins.push_back(19.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));

  ExcitEnergies.push_back(11825*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(70*keV));
}

G4C15GEMProbability::~G4C15GEMProbability()
{}
