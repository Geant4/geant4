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
// $Id: G4Na23GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Na23GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Na23GEMProbability::G4Na23GEMProbability() :
  G4GEMProbability(23,11,3.0/2.0) // A,Z,Spin
{

    ExcitEnergies.push_back(439.80*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(1109.0e-3*picosecond);

    ExcitEnergies.push_back(2076.4*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(28.0e-3*picosecond);

    ExcitEnergies.push_back(2390.9*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(555.0e-3*picosecond);

    ExcitEnergies.push_back(2639.8*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(69.0e-3*picosecond);

    ExcitEnergies.push_back(2703.7*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(76.0e-3*picosecond);

    ExcitEnergies.push_back(2982.4*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(3.5e-3*picosecond);

    ExcitEnergies.push_back(3678.3*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(22.0e-3*picosecond);

    ExcitEnergies.push_back(3848.2*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(90.0e-3*picosecond);

    ExcitEnergies.push_back(3914.7*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(7.6e-3*picosecond);

    ExcitEnergies.push_back(4432.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(0.19e-3*picosecond);

    ExcitEnergies.push_back(4775.6*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(1.4e-3*picosecond);

    ExcitEnergies.push_back(5377.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.25e-3*picosecond);

    ExcitEnergies.push_back(5533.0*keV);
    ExcitSpins.push_back(11.0/2.0);
    ExcitLifetimes.push_back(8.0e-3*picosecond);

    ExcitEnergies.push_back(6117.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(52.0e-3*picosecond);

    ExcitEnergies.push_back(6237.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(6348.0*keV);
    ExcitSpins.push_back(9.0/2.0);
    ExcitLifetimes.push_back(26e-3*picosecond);

    ExcitEnergies.push_back(7080.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(0.24e-3*picosecond);

    ExcitEnergies.push_back(7132.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.21e-3*picosecond);

    ExcitEnergies.push_back(7888.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(0.17e-3*picosecond);

    ExcitEnergies.push_back(8663.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(0.10e-3*picosecond);

    ExcitEnergies.push_back(8945.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(21.0e-3*picosecond);

    ExcitEnergies.push_back(9209.6*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.11*keV));
    
    ExcitEnergies.push_back(9250.8*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1.6*keV));
    
    ExcitEnergies.push_back(9395.4*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.03*keV));
    
    ExcitEnergies.push_back(9404.4*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(65.0*keV));
    
    ExcitEnergies.push_back(9425.6*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.6*keV));
    
    ExcitEnergies.push_back(9485.9*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.22*keV));
    
    ExcitEnergies.push_back(9606.6*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(6.0*keV));
    
    ExcitEnergies.push_back(9650.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(2.1*keV));

    ExcitEnergies.push_back(9653.7*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(105.0*keV));
    
    ExcitEnergies.push_back(9679.5*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.6*keV));
    
    ExcitEnergies.push_back(9698.9*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(29.0*keV));
    
    ExcitEnergies.push_back(9730.8*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(0.35*keV));
    
    ExcitEnergies.push_back(9753.5*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(4.0*keV));
    
    ExcitEnergies.push_back(9813.8*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1.0*keV));
    
    ExcitEnergies.push_back(9833.4*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(47.0*keV));
    
    ExcitEnergies.push_back(9849.1*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(150.0*keV));

    ExcitEnergies.push_back(9998.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(475.0*keV));
    
    ExcitEnergies.push_back(10014.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(69.0*keV));
    
    ExcitEnergies.push_back(10083.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1270.0*keV));
    
    ExcitEnergies.push_back(10110.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(4200.0*keV));
    
    ExcitEnergies.push_back(10167.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(65.0*keV));
    
    ExcitEnergies.push_back(10229.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(5.0*keV));
    
    ExcitEnergies.push_back(10239.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(2450.0*keV));
    
    ExcitEnergies.push_back(10313.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(2000.0*keV));
    
    ExcitEnergies.push_back(10334.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(190.0*keV));
    
    ExcitEnergies.push_back(10341.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(11.0*keV));

    ExcitEnergies.push_back(10348.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(210.0*keV));
    
    ExcitEnergies.push_back(10436.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(25.0*keV));
    
    ExcitEnergies.push_back(10476.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(470.0*keV));
    
    ExcitEnergies.push_back(10499.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(920.0*keV));
    
    ExcitEnergies.push_back(10505.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(560.0*keV));
    
    ExcitEnergies.push_back(10516.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(100.0*keV));
    
    ExcitEnergies.push_back(10546.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(540.0*keV));
    
    ExcitEnergies.push_back(10570.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1100.0*keV));
    
    ExcitEnergies.push_back(10615.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(425.0*keV));
    
    ExcitEnergies.push_back(10673.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(23000.0*keV));

    ExcitEnergies.push_back(10702.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(400.0*keV));
    
    ExcitEnergies.push_back(10769.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(5.0*keV));
    
    ExcitEnergies.push_back(10822.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(1700.0*keV));
    
    ExcitEnergies.push_back(10823.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(26000.0*keV));
    
    ExcitEnergies.push_back(10837.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(100.0*keV));
    
    ExcitEnergies.push_back(10869.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(21000.0*keV));
    
    ExcitEnergies.push_back(10904.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(53.0*keV));
    
    ExcitEnergies.push_back(10907.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(2850.0*keV));
    
    ExcitEnergies.push_back(10908.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(900.0*keV));
    
    ExcitEnergies.push_back(10917.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(55.0*keV));

    ExcitEnergies.push_back(10931.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(3500.0*keV));
    
    ExcitEnergies.push_back(10949.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(5200.0*keV));
    
    ExcitEnergies.push_back(10954.0*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(65.0*keV));
    
    ExcitEnergies.push_back(10972.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(400.0*keV));
    
    ExcitEnergies.push_back(10975.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(18.0*keV));
    
    ExcitEnergies.push_back(10979.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(6000.0*keV));
    
    ExcitEnergies.push_back(10993.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(20600.0*keV));
    
    ExcitEnergies.push_back(10994.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(60.0*keV));
    
    ExcitEnergies.push_back(11039.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(500.0*keV));
    
    ExcitEnergies.push_back(11087.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(800.0*keV));
    
    ExcitEnergies.push_back(11109.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(135.0*keV));
    
    ExcitEnergies.push_back(11110.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(fPlanck/(4100.0*keV));
    
}

G4Na23GEMProbability::~G4Na23GEMProbability() 
{}
