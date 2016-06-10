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
//
// $Id: G4Be10GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Be10GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Be10GEMProbability::G4Be10GEMProbability() :
  G4GEMProbability(10,4,0.0) // A,Z,Spin
{
  ExcitEnergies.push_back(3368.03*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(125.0e-15*s);

  ExcitEnergies.push_back(5958.39*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(55.0e-15*s);

  ExcitEnergies.push_back(5959.9*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(0.8e-12*s);

  ExcitEnergies.push_back(7371.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(15.7*keV));

  ExcitEnergies.push_back(7542.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(6.3*keV));

  ExcitEnergies.push_back(9270.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(9400.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(291.0*keV));

  ExcitEnergies.push_back(11.76e3*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(121.0*keV));

  ExcitEnergies.push_back(17790.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(110.0*keV));

  ExcitEnergies.push_back(18550.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(350.0*keV));
}

G4Be10GEMProbability::~G4Be10GEMProbability()
{}
