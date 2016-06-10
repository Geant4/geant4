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
// $Id: G4Ne21GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Ne21GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Ne21GEMProbability::G4Ne21GEMProbability() :
  G4GEMProbability(21,10,3.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(350.72*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(7.09*picosecond);

  ExcitEnergies.push_back(1745.6*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(55.0e-3*picosecond);

  ExcitEnergies.push_back(2788.5*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(81.0*picosecond);

  ExcitEnergies.push_back(2796.1*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(7.0e-3*picosecond);

  ExcitEnergies.push_back(2865.6*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(35.0e-3*picosecond);

  ExcitEnergies.push_back(3662.1*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(59.0e-3*picosecond);

  ExcitEnergies.push_back(3733.7*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(14.0e-3*picosecond);

  ExcitEnergies.push_back(3882.9*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(30.0e-3*picosecond);

  ExcitEnergies.push_back(4432.2*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(25.0e-3*picosecond);

  ExcitEnergies.push_back(4524.2*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(7.0e-3*picosecond);

  ExcitEnergies.push_back(4683.6*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(11.0e-3*picosecond);

  ExcitEnergies.push_back(4725.7*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(7.0e-3*picosecond);

  ExcitEnergies.push_back(5334.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(7.0e-3*picosecond);

  ExcitEnergies.push_back(5430.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(14.0e-3*picosecond);

  ExcitEnergies.push_back(5525.0*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(69.0e-3*picosecond);

  ExcitEnergies.push_back(5550.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(28.0e-3*picosecond);

  ExcitEnergies.push_back(5629.4*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(7.0e-3*picosecond);

  ExcitEnergies.push_back(5690.5*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(7.0e-3*picosecond);

  ExcitEnergies.push_back(5775.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(28.0e-3*picosecond);

  ExcitEnergies.push_back(5821.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(55.0e-3*picosecond);

  ExcitEnergies.push_back(5823.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(24.0e-3*picosecond);

  ExcitEnergies.push_back(5992.9*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(7.0e-3*picosecond);

  ExcitEnergies.push_back(6030.7*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(24.0e-7*picosecond);

  ExcitEnergies.push_back(6169.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(24.0e-3*picosecond);

  ExcitEnergies.push_back(6265.1*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(24.0e-3*picosecond);

  ExcitEnergies.push_back(6446.6*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(14.0e-3*picosecond);

  ExcitEnergies.push_back(6553.0*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(31.0e-3*picosecond);

  ExcitEnergies.push_back(6605.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(24.0e-3*picosecond);

  ExcitEnergies.push_back(6642.0*keV);
  ExcitSpins.push_back(11.0/2.0);
  ExcitLifetimes.push_back(66.0e-3*picosecond);

  ExcitEnergies.push_back(6747.4*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(28.0e-3*picosecond);

  ExcitEnergies.push_back(7212.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(107.0*keV));

  ExcitEnergies.push_back(7653.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(14.0*keV));

  ExcitEnergies.push_back(7980.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6.0*keV));

  ExcitEnergies.push_back(8008.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(32.0*keV));

  ExcitEnergies.push_back(8062.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(8.0*keV));

  ExcitEnergies.push_back(8281.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(27.0*keV));

  ExcitEnergies.push_back(8352.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(10.0*keV));

  ExcitEnergies.push_back(8583.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(38.0*keV));

  ExcitEnergies.push_back(8660.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(54.0*keV));

  ExcitEnergies.push_back(8781.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(50*keV));

  ExcitEnergies.push_back(8857.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2.8*keV));

}

G4Ne21GEMProbability::~G4Ne21GEMProbability() 
{}

