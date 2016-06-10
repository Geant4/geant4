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
// $Id: G4Mg25GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Mg25GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Mg25GEMProbability::G4Mg25GEMProbability() :
  G4GEMProbability(25,12,5.0/2.0) // A,Z,Spin
{

    ExcitEnergies.push_back(585.09*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(3.38*nanosecond);

    ExcitEnergies.push_back(974.79*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(11.4*picosecond);

    ExcitEnergies.push_back(1611.8*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(14.6e-3*picosecond);

    ExcitEnergies.push_back(1964.7*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.69*picosecond);

    ExcitEnergies.push_back(2563.8*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(2737.7*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(277.0e-3*picosecond);

    ExcitEnergies.push_back(2801.1*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(28.0e-3*picosecond);

    ExcitEnergies.push_back(3405.2*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(7.1e-3*picosecond);

    ExcitEnergies.push_back(3413.7*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(3907.8*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(3970.7*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(19.0e-3*picosecond);

    ExcitEnergies.push_back(4059.6*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(53.0e-3*picosecond);

    ExcitEnergies.push_back(4277.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(4359.4*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(4711.4*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(26.0e-3*picosecond);

    ExcitEnergies.push_back(4722.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(5012.2*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(5251.3*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(5461.7*keV);
    ExcitSpins.push_back(13.0/2.0);
    ExcitLifetimes.push_back(1.7*picosecond);

    ExcitEnergies.push_back(5474.8*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(5520.9*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(5533.6*keV);
    ExcitSpins.push_back(11.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(5746.6*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(0.58e-3*picosecond);

    ExcitEnergies.push_back(5793.2*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(45.0e-3*picosecond);

    ExcitEnergies.push_back(5859.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(5971.6*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(5978.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(6041.2*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(6082.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(6775.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(7088.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(7282.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(7224.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(6.9e-3*picosecond);

    ExcitEnergies.push_back(7411.9*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(7.6*keV));

    ExcitEnergies.push_back(7587.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(78.0*keV));

    ExcitEnergies.push_back(7744.9*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(31.0*keV));

    ExcitEnergies.push_back(7787.9*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(14.0*eV));

    ExcitEnergies.push_back(7809.8*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.53*keV));

    ExcitEnergies.push_back(7864.5*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(14.0*eV));

    ExcitEnergies.push_back(7947.7*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.9*keV));

    ExcitEnergies.push_back(7964.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(21.0*keV));

    ExcitEnergies.push_back(8141.8*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(13.0*keV));

    ExcitEnergies.push_back(8312.3*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1.2*keV));

    ExcitEnergies.push_back(8325.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(10.0*keV));

    ExcitEnergies.push_back(8363.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(20.0*keV));

    ExcitEnergies.push_back(8559.8*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1.8*keV));

    ExcitEnergies.push_back(8570.5*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(4.2*keV));

    ExcitEnergies.push_back(8580.5*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(8.4*keV));

    ExcitEnergies.push_back(8835.5*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(4.1*keV));

    ExcitEnergies.push_back(8870.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(15.0*keV));

    ExcitEnergies.push_back(8895.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(35.0*keV));

    ExcitEnergies.push_back(8972.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(9.6*keV));
}

G4Mg25GEMProbability::~G4Mg25GEMProbability() 
{}


