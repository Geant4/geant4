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
// $Id: G4C13GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4C13GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4C13GEMProbability::G4C13GEMProbability() :
  G4GEMProbability(13,6,1.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(3089.443*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(1.04e-15*s);

  ExcitEnergies.push_back(3684.507*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(1.04e-15*s);

  ExcitEnergies.push_back( 3853.807*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(7.5e-12*s);

  ExcitEnergies.push_back(6864*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6*keV));

  ExcitEnergies.push_back(7492*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5*keV));

  ExcitEnergies.push_back(7547*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.2*keV));

  ExcitEnergies.push_back(7677*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(70*keV));

  ExcitEnergies.push_back(8.2E3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000*keV));

  ExcitEnergies.push_back(8860*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(150*keV));

  ExcitEnergies.push_back(9498*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5*keV));

  ExcitEnergies.push_back(9897*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(26*keV));

  ExcitEnergies.push_back(10753*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(55*keV));

  ExcitEnergies.push_back(10818*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(24*keV));

  ExcitEnergies.push_back(10996*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(37*keV));

  ExcitEnergies.push_back(11080*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4*keV));

  ExcitEnergies.push_back(11851*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(68*keV));

  ExcitEnergies.push_back(11970*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(200*keV));

  ExcitEnergies.push_back(12106*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(81*keV));

  ExcitEnergies.push_back(12400*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(150*keV));

  ExcitEnergies.push_back(13280*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(340*keV));

  ExcitEnergies.push_back(13410*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(35*keV));

  ExcitEnergies.push_back(13560*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(500*keV));

  ExcitEnergies.push_back(13760*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(300*keV));

  ExcitEnergies.push_back(14120*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(200*keV));

  ExcitEnergies.push_back(14.39E3*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(260*keV));

  ExcitEnergies.push_back(14940*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(380*keV));

  ExcitEnergies.push_back(15108.2*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5*keV));

  ExcitEnergies.push_back(19500*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(450*keV));

}

G4C13GEMProbability::~G4C13GEMProbability() 
{}
