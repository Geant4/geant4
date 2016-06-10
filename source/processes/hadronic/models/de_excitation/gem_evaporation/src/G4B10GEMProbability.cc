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
// $Id: G4B10GEMProbability.cc 86783 2014-11-18 08:43:58Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4B10GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4B10GEMProbability::G4B10GEMProbability() :
  G4GEMProbability(10,5,3.0) // A,Z,Spin
{
    ExcitEnergies.push_back(718.38*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(0.707e-9*s);

    ExcitEnergies.push_back(1740.05*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(5.0e-15*s);

    ExcitEnergies.push_back(2154.27*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(1.48e-12*s);

    ExcitEnergies.push_back(3587.13*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(106.0e-15*s);

    ExcitEnergies.push_back(4774.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(8.7*keV));

    ExcitEnergies.push_back(5110.3*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(0.98*keV));

    ExcitEnergies.push_back(5163.9*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(4.2e-15*s);

    ExcitEnergies.push_back(5182.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(110.0*keV));

    ExcitEnergies.push_back(5919.5*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(6.0*keV));

    ExcitEnergies.push_back(6024.9*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(fPlanck/(0.05*keV));

    ExcitEnergies.push_back(6127.2*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(2.36*keV));

    ExcitEnergies.push_back(6561.0*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(fPlanck/(25.1*keV));

    ExcitEnergies.push_back(6875.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(120.0*keV));

    ExcitEnergies.push_back(7002.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(100.0*keV));

    ExcitEnergies.push_back(7428.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(100.0*keV));

    ExcitEnergies.push_back(7467.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(65.0*keV));

    ExcitEnergies.push_back(7479.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(74.0*keV));

    ExcitEnergies.push_back(7559.9*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(2.65*keV));

    ExcitEnergies.push_back(7750.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(250.0*keV));

    ExcitEnergies.push_back(7819.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(260.0*keV));

    ExcitEnergies.push_back(8070.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(800.0*keV));

    ExcitEnergies.push_back(8700.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(200.0*keV));

    ExcitEnergies.push_back(8889.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(84.0*keV));

    ExcitEnergies.push_back(8895.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(40.0*keV));

    ExcitEnergies.push_back(9700.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(700.0*keV));

    ExcitEnergies.push_back(10840.0*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(fPlanck/(300.0*keV));

    ExcitEnergies.push_back(11520.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(500.0*keV));

    ExcitEnergies.push_back(12560.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(100.0*keV));

    ExcitEnergies.push_back(13490.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(300.0*keV));

    ExcitEnergies.push_back(14.4e3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(800.0*keV));

    ExcitEnergies.push_back(18.2e3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(1500.0*keV));

    ExcitEnergies.push_back(18430.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(340.0*keV));

    ExcitEnergies.push_back(18800.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(600.0*keV));

    ExcitEnergies.push_back(19290.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(190.0*keV));
}

G4B10GEMProbability::~G4B10GEMProbability()
{}
