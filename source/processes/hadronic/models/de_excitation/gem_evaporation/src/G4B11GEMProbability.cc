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
// $Id: G4B11GEMProbability.cc 86783 2014-11-18 08:43:58Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4B11GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4B11GEMProbability::G4B11GEMProbability() :
  G4GEMProbability(11,5,3.0/2.0) // A,Z,Spin
{
    ExcitEnergies.push_back(2124.693*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(3.8e-15*s);

    ExcitEnergies.push_back(4444.98*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.82e-15*s);

    ExcitEnergies.push_back(5020.3*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(0.236e-15*s);

    ExcitEnergies.push_back(6741.85*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(15.0e-15*s);

    ExcitEnergies.push_back(6791.8*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(1.18e-15*s);

    ExcitEnergies.push_back(7285.51*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.4e-15*s);

    ExcitEnergies.push_back(7977.84*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(0.4e-15*s);

    ExcitEnergies.push_back(8560.1*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.49e-15*s);

    ExcitEnergies.push_back(8920.47*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(4.4e-15*s);

    ExcitEnergies.push_back(9183.5*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1.9*eV));

    ExcitEnergies.push_back(9271.7*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(4.0*keV));

    ExcitEnergies.push_back(9876.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(110.0*keV));

    ExcitEnergies.push_back(10260.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(165.0*keV));

    ExcitEnergies.push_back(10330.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(110.0*keV));

    ExcitEnergies.push_back(10597.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(100.0*keV));

    ExcitEnergies.push_back(10960.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(4500.0*keV));

    ExcitEnergies.push_back(11265.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(110.0*keV));

    ExcitEnergies.push_back(11444.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(103.0*keV));

    ExcitEnergies.push_back(11886.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(200.0*keV));

    ExcitEnergies.push_back(12.0e3*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

    ExcitEnergies.push_back(12557.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(210.0*keV));

    ExcitEnergies.push_back(12916.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(155.0*keV));

    ExcitEnergies.push_back(13137.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(426.0*keV));

    ExcitEnergies.push_back(13.16e3*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(430.0*keV));

    ExcitEnergies.push_back(14.04e3*keV);
    ExcitSpins.push_back(11.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(500.0*keV));

    ExcitEnergies.push_back(14.34e3*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(254.0*keV));

    ExcitEnergies.push_back(14565.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(30.0*keV));

    ExcitEnergies.push_back(15.32e3*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(635.0*keV));

    ExcitEnergies.push_back(16437.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(30.0*keV));
}

G4B11GEMProbability::~G4B11GEMProbability()
{}
