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
// $Id: G4F19GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4F19GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4F19GEMProbability::G4F19GEMProbability() :
  G4GEMProbability(19,9,1.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(109.894*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(0.580e-3*picosecond);

  ExcitEnergies.push_back(197.143*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(89.30e-3*picosecond);

  ExcitEnergies.push_back(1345.67*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(3.3*picosecond);

  ExcitEnergies.push_back(1458.7*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(54.0e-3*picosecond);

  ExcitEnergies.push_back(1554.038*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(3.0e-3*picosecond);

  ExcitEnergies.push_back(2779.849*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(181.0e-3*picosecond);

  ExcitEnergies.push_back(3908.17*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(14.0e-3*picosecond);

  ExcitEnergies.push_back(3998.7*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(16.0e-3*picosecond);

  ExcitEnergies.push_back(4032.5*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(51.0e-3*picosecond);

  ExcitEnergies.push_back(4555.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(14.0e-3*picosecond);

  ExcitEnergies.push_back(4557.5*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(14.0e-3*picosecond);

  ExcitEnergies.push_back(4648.0*keV);
  ExcitSpins.push_back(13.0/2.0);
  ExcitLifetimes.push_back(1.5*picosecond);

  ExcitEnergies.push_back(4683.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(10.7e-3*picosecond);

  ExcitEnergies.push_back(5340.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(10.0e-3*picosecond);

  ExcitEnergies.push_back(5464.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(13.0e-3*picosecond);

  ExcitEnergies.push_back(6076.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.2*keV));

  ExcitEnergies.push_back(6093.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(6250.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(8.0*keV));

  ExcitEnergies.push_back(6290.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2.4*keV));

  ExcitEnergies.push_back(6332.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2.4*keV));

  ExcitEnergies.push_back(6430.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(280.0*keV));

  ExcitEnergies.push_back(6525.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(6553.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.6*keV));

  ExcitEnergies.push_back(6788.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2.4*keV));

  ExcitEnergies.push_back(6838.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(1.2*keV));

  ExcitEnergies.push_back(6890.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(28.0*keV));

  ExcitEnergies.push_back(6926.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2.4*keV));

  ExcitEnergies.push_back(6990.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(51.0*keV));

  ExcitEnergies.push_back(7110.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(32.0*keV));

  ExcitEnergies.push_back(7120.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(8.0*keV));

  ExcitEnergies.push_back(7364.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(63.0*keV));

  ExcitEnergies.push_back(7702.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(24.0*keV));

  ExcitEnergies.push_back(8014.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6.0*keV));

  ExcitEnergies.push_back(8086.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6.0*keV));

  ExcitEnergies.push_back(8136.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5.0*keV));

  ExcitEnergies.push_back(8195.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(8.0*keV));

  ExcitEnergies.push_back(8590.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(2.0*keV));

  ExcitEnergies.push_back(8637.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(95.0*keV));

  ExcitEnergies.push_back(8795.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(45.0*keV));

  ExcitEnergies.push_back(8928.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.6*keV));

  ExcitEnergies.push_back(9098.4*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(20.0e-3*keV));

  ExcitEnergies.push_back(9166.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5.8*keV));

  ExcitEnergies.push_back(9321.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.9*keV));

  ExcitEnergies.push_back(9527.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(29.0*keV));

  ExcitEnergies.push_back(9578.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(26.0*keV));

  ExcitEnergies.push_back(9668.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(3.8*keV));

  ExcitEnergies.push_back(9819.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(0.29*keV));

  ExcitEnergies.push_back(9888.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(31.0*keV));

  ExcitEnergies.push_back(10136.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.7*keV));

  ExcitEnergies.push_back(10161.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(31.0*keV));

  ExcitEnergies.push_back(10231.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.3*keV));

  ExcitEnergies.push_back(10253.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(23.0*keV));

  ExcitEnergies.push_back(10306.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(9.2*keV));

  ExcitEnergies.push_back(10496.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(6.2*keV));

  ExcitEnergies.push_back(10554.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(7.6*keV));

  ExcitEnergies.push_back(10580.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(22.0*keV));

  ExcitEnergies.push_back(10613.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(4.8*keV));

  ExcitEnergies.push_back(10763.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(5.4*keV));

  ExcitEnergies.push_back(10858.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(25.0*keV));

  ExcitEnergies.push_back(10972.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(11.0*keV));

  ExcitEnergies.push_back(11070.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(31.0*keV));

  ExcitEnergies.push_back(11199.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(fPlanck/(43.0*keV));

}

G4F19GEMProbability::~G4F19GEMProbability() 
{}
