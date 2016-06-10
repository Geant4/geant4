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
// $Id: G4Ne18GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Ne18GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Ne18GEMProbability::G4Ne18GEMProbability() :
  G4GEMProbability(18,10,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1887.3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.34*picosecond);

  ExcitEnergies.push_back(3376.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(3.0*picosecond);

  ExcitEnergies.push_back(3576.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(1.4*picosecond);

  ExcitEnergies.push_back(3616.4*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.040*picosecond);

  ExcitEnergies.push_back(4510.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(4580.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(7062.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(180.0*keV));

  ExcitEnergies.push_back(7915.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

}

G4Ne18GEMProbability::~G4Ne18GEMProbability()
{}
