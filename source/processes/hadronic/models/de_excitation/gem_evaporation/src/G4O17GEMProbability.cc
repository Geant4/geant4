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
// $Id: G4O17GEMProbability.cc 74869 2013-10-23 09:26:17Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4O17GEMProbability.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4Log.hh"

G4O17GEMProbability::G4O17GEMProbability() :
  G4GEMProbability(17,9,5.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(870.81*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(179.2*picosecond);

  ExcitEnergies.push_back(3055.2*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(83.0E-3*picosecond);

  ExcitEnergies.push_back(3841.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(17.0E-3*picosecond);

  ExcitEnergies.push_back(4553.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(5086.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(95.0*keV));

  ExcitEnergies.push_back(5380.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(28.0*keV));

  ExcitEnergies.push_back(5698.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.4*keV));

  ExcitEnergies.push_back(5870.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6.6*keV));

  ExcitEnergies.push_back(5940.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(32.0*keV));

  ExcitEnergies.push_back(6357.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(124.0*keV));

  ExcitEnergies.push_back(7168.7*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.5*keV));

  ExcitEnergies.push_back(7202.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(280.0*keV));

  ExcitEnergies.push_back(7383.1*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(0.6*keV));

  ExcitEnergies.push_back(7386.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(0.9*keV));

  ExcitEnergies.push_back(7560.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(500.0*keV));

  ExcitEnergies.push_back(7577.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.0*keV));

  ExcitEnergies.push_back(7690.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(18.0*keV));

  ExcitEnergies.push_back(7956.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(90.0*keV));

  ExcitEnergies.push_back(7990.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(270.0*keV));

  ExcitEnergies.push_back(8070.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(85.0*keV));

  ExcitEnergies.push_back(8180.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(69.0*keV));

  ExcitEnergies.push_back(8200.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(8352.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(9.0*keV));

  ExcitEnergies.push_back(8410.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(8474.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(7.0*keV));

  ExcitEnergies.push_back(8508.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*G4Log(2.0)/(5.0*keV));

  ExcitEnergies.push_back(8700.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

  ExcitEnergies.push_back(8898.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(101.0*keV));

  ExcitEnergies.push_back(8972.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(21.0*keV));

  ExcitEnergies.push_back(9148.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(9187.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.0*keV));

  ExcitEnergies.push_back(9201.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5.5*keV));

  ExcitEnergies.push_back(9422.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(120.0*keV));

  ExcitEnergies.push_back(9493.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*G4Log(2.0)/(15.0*keV));

  ExcitEnergies.push_back(9720.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(16.0*keV));

  ExcitEnergies.push_back(9775.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(25.0*keV));

  ExcitEnergies.push_back(9977.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(10178.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(10337.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(10490.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(75.0*keV));

  ExcitEnergies.push_back(10773.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(10910.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(57.0*keV));

  ExcitEnergies.push_back(11076.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5.0*keV));

  ExcitEnergies.push_back(12946.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6.0*keV));

  ExcitEnergies.push_back(14980.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(21.7E3*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(750.0*keV));

  ExcitEnergies.push_back(22.1E3*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(750.0*keV));

  ExcitEnergies.push_back(22.5E3*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

  ExcitEnergies.push_back(23.0E3*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(400.0*keV));
}

G4O17GEMProbability::~G4O17GEMProbability() 
{}
