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
// $Id: G4C10GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4C10GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4C10GEMProbability::G4C10GEMProbability() :
  G4GEMProbability(10,6,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(3353.7*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(107.0e-15*s);

  ExcitEnergies.push_back(6580.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));
}

G4C10GEMProbability::~G4C10GEMProbability()
{}
