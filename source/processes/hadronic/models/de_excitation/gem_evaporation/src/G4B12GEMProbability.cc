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
// $Id: G4B12GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4B12GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4B12GEMProbability::G4B12GEMProbability() :
  G4GEMProbability(12,5,1.0) // A,Z,Spin
{
    ExcitEnergies.push_back(953.14*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(180.0e-15*s);

    ExcitEnergies.push_back(1673.65*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(35.0e-15*s);

    ExcitEnergies.push_back(2620.8*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(49.0e-15*s);

    ExcitEnergies.push_back(3388.3*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(3.1*eV));

    ExcitEnergies.push_back(3759.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(40*keV));

    ExcitEnergies.push_back(4301*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(9*keV));

    ExcitEnergies.push_back(4518*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(fPlanck/(110*keV));

    ExcitEnergies.push_back(5.00E+3*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(50*keV));

    ExcitEnergies.push_back(5612*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(110*keV));

    ExcitEnergies.push_back(5726*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(50*keV));

    ExcitEnergies.push_back(6.6E+3*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(140*keV));

    ExcitEnergies.push_back(7.67E+3*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(45*keV));

    ExcitEnergies.push_back(7836*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(60*keV));

    ExcitEnergies.push_back(7937*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(27*keV));

    ExcitEnergies.push_back(8.24E+3*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(65*keV));

    ExcitEnergies.push_back(8.58E+3*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(75*keV));

    ExcitEnergies.push_back(9040*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(fPlanck/(95*keV));

    ExcitEnergies.push_back(9585*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(fPlanck/(34*keV));

    ExcitEnergies.push_back(1.275E+4*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(fPlanck/(85*keV));

    ExcitEnergies.push_back(1.482E+4*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(fPlanck/(200*keV));
}

G4B12GEMProbability::~G4B12GEMProbability()
{}

