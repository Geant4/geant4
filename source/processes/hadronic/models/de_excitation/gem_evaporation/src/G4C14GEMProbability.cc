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
// $Id: G4C14GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4C14GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4C14GEMProbability::G4C14GEMProbability() :
  G4GEMProbability(14,6,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(6093.8*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(7e-15*s);

  ExcitEnergies.push_back(6589.4*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(0.4e-12*s);

  ExcitEnergies.push_back(6728.2*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(67e-12*s);

  ExcitEnergies.push_back(6902.6*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(25e-15*s);

  ExcitEnergies.push_back(7012.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(9e-15*s);

  ExcitEnergies.push_back(8318.3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(3.4*keV));

  ExcitEnergies.push_back(9799*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(45*keV));

  ExcitEnergies.push_back(10437*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(16*keV));

  ExcitEnergies.push_back(10509*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(26*keV));

  ExcitEnergies.push_back(11306*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(46*keV));

  ExcitEnergies.push_back(11397*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(22*keV));

  ExcitEnergies.push_back(11667*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(20*keV));

  ExcitEnergies.push_back(12860*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));

  ExcitEnergies.push_back(12964*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(30*keV));

  ExcitEnergies.push_back(14667*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(57*keV));
}

G4C14GEMProbability::~G4C14GEMProbability()
{}
