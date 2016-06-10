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
// $Id: G4Ne20GEMProbability.cc 67983 2013-03-13 10:42:03Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4Ne20GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4Ne20GEMProbability::G4Ne20GEMProbability() :
  G4GEMProbability(20,10,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1633.8*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.83*picosecond);

  ExcitEnergies.push_back(4247.3*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(64.0e-3*picosecond);

  ExcitEnergies.push_back(4968.2*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(3.3*picosecond);

  ExcitEnergies.push_back(5621.7*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(139.0e-3*picosecond);

  ExcitEnergies.push_back(5785.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(1.3E-2*keV));

  ExcitEnergies.push_back(6722.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(15.0*keV));


  ExcitEnergies.push_back(7005.5*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(305.0e-3*picosecond);


  ExcitEnergies.push_back(7166.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(8.0*keV));

  ExcitEnergies.push_back(7196.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(7424.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(8.0*keV));

  ExcitEnergies.push_back(7834.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(2.0*keV));

  ExcitEnergies.push_back(8600.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(800.0*keV));

  ExcitEnergies.push_back(8720.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(2.5*keV));

  ExcitEnergies.push_back(8775.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(0.110*keV));

  ExcitEnergies.push_back(8800.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(800.0*keV));

  ExcitEnergies.push_back(8820.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(1.0*keV));

  ExcitEnergies.push_back(8850.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(19.0*keV));

  ExcitEnergies.push_back(9040.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(3.0*keV));

  ExcitEnergies.push_back(9117.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(3.2*keV));

  ExcitEnergies.push_back(9489.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(29.0*keV));

  ExcitEnergies.push_back(9950.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(24.0*keV));

  ExcitEnergies.push_back(9990.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(10257.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(141.0*keV));

  ExcitEnergies.push_back(10260.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(2.0*keV));

  ExcitEnergies.push_back(10401.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(81.0*keV));

  ExcitEnergies.push_back(10548.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(16.0*keV));

  ExcitEnergies.push_back(10579.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(24.0*keV));

  ExcitEnergies.push_back(10609.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(16.0e-3*picosecond);

  ExcitEnergies.push_back(10790.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(350.0*keV));

  ExcitEnergies.push_back(10836.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(13.0*keV));

  ExcitEnergies.push_back(10836.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(45.0*keV));

  ExcitEnergies.push_back(10970.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(580.0*keV));

  ExcitEnergies.push_back(11015.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(24.0*keV));

  ExcitEnergies.push_back(11080.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(3.0*keV));

  ExcitEnergies.push_back(11230.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(172.0*keV));

  ExcitEnergies.push_back(11270.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(11324.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(53.0*keV));

  ExcitEnergies.push_back(11871.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(46.0*keV));

  ExcitEnergies.push_back(11925.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(0.44*keV));

  ExcitEnergies.push_back(11948.0*keV);
  ExcitSpins.push_back(8.0);
  ExcitLifetimes.push_back(fPlanck/(35.0e-3*keV));

  ExcitEnergies.push_back(11953.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(24.0*keV));

  ExcitEnergies.push_back(11971.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(29.0*keV));

  ExcitEnergies.push_back(12150.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(12224.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(142.0*keV));

  ExcitEnergies.push_back(12245.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(12367.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(46.0*keV));

  ExcitEnergies.push_back(12410.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(8.0*keV));

  ExcitEnergies.push_back(12559.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(101.0*keV));

  ExcitEnergies.push_back(12770.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(12980.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(13060.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(1.0*keV));

  ExcitEnergies.push_back(13086.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(13168.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(2.3*keV));

  ExcitEnergies.push_back(13180.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(13224.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(95.0*keV));

  ExcitEnergies.push_back(13224.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(95.0*keV));

  ExcitEnergies.push_back(13304.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(0.9*keV));

  ExcitEnergies.push_back(13333.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(80.0e-3*keV));

  ExcitEnergies.push_back(13342.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(20.0*keV));

  ExcitEnergies.push_back(13411.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(35.0*keV));

  ExcitEnergies.push_back(13420.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(110.0*keV));

  ExcitEnergies.push_back(13462.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(190.0*keV));

  ExcitEnergies.push_back(13479.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(7.1*keV));

  ExcitEnergies.push_back(13523.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(30.0*keV));

  ExcitEnergies.push_back(13541.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(63.0*keV));

  ExcitEnergies.push_back(13584.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(10.0*keV));

  ExcitEnergies.push_back(13650.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(22.0*keV));

  ExcitEnergies.push_back(13660.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(115.0*keV));

  ExcitEnergies.push_back(13673.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(4.5*keV));

  ExcitEnergies.push_back(13700.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(320.0*keV));

  ExcitEnergies.push_back(13730.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(170.0*keV));

  ExcitEnergies.push_back(13733.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(7.7*keV));

  ExcitEnergies.push_back(13870.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(190.0*keV));

  ExcitEnergies.push_back(13880.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(13903.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(48.0*keV));

  ExcitEnergies.push_back(13946.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(14017.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(14030.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(140.0*keV));

  ExcitEnergies.push_back(14124.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(4.7*keV));

  ExcitEnergies.push_back(14134.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(51.0*keV));

  ExcitEnergies.push_back(14148.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(11.8*keV));

  ExcitEnergies.push_back(14197.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(13.9*keV));

  ExcitEnergies.push_back(14300.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(240.0*keV));

  ExcitEnergies.push_back(14467.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(135.0*keV));

  ExcitEnergies.push_back(14600.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(240.0*keV));

  ExcitEnergies.push_back(14604.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(125.0*keV));

  ExcitEnergies.push_back(14695.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(38.0*keV));

  ExcitEnergies.push_back(14850.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(15030.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(90.0*keV));

  ExcitEnergies.push_back(15260.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(285.0*keV));

  ExcitEnergies.push_back(15300.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(285.0*keV));

  ExcitEnergies.push_back(15618.0*keV);
  ExcitSpins.push_back(8.0);
  ExcitLifetimes.push_back(fPlanck/(28.0*keV));

  ExcitEnergies.push_back(16728.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(fPlanck/(2.0*keV));

  ExcitEnergies.push_back(18080.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(140.0*keV));

  ExcitEnergies.push_back(18310.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(240.0*keV));

  ExcitEnergies.push_back(18426.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(10.0*keV));

  ExcitEnergies.push_back(18700.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(600.0*keV));

  ExcitEnergies.push_back(19160.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(19400.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(280.0*keV));

  ExcitEnergies.push_back(19840.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(280.0*keV));

  ExcitEnergies.push_back(20160.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(250.0*keV));

  ExcitEnergies.push_back(20400.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(360.0*keV));

  ExcitEnergies.push_back(20400.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(20680.0*keV);
  ExcitSpins.push_back(9.0);
  ExcitLifetimes.push_back(fPlanck/(120.0*keV));

  ExcitEnergies.push_back(21000.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(21080.0*keV);
  ExcitSpins.push_back(9.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(21300.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(300.0*keV));

  ExcitEnergies.push_back(21800.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(300.0*keV));

  ExcitEnergies.push_back(22300.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(500.0*keV));

  ExcitEnergies.push_back(22700.0*keV);
  ExcitSpins.push_back(9.0);
  ExcitLifetimes.push_back(fPlanck/(500.0*keV));

  ExcitEnergies.push_back(22840.0*keV);
  ExcitSpins.push_back(9.0);
  ExcitLifetimes.push_back(fPlanck/(250.0*keV));

  ExcitEnergies.push_back(23400.0*keV);
  ExcitSpins.push_back(8.0);
  ExcitLifetimes.push_back(fPlanck/(500.0*keV));

  ExcitEnergies.push_back(24110.0*keV);
  ExcitSpins.push_back(8.0);
  ExcitLifetimes.push_back(fPlanck/(350.0*keV));

  ExcitEnergies.push_back(25000.0*keV);
  ExcitSpins.push_back(8.0);
  ExcitLifetimes.push_back(fPlanck/(600.0*keV));

  ExcitEnergies.push_back(28000.0*keV);
  ExcitSpins.push_back(8.0);
  ExcitLifetimes.push_back(fPlanck/(1600.0*keV));

}

G4Ne20GEMProbability::~G4Ne20GEMProbability() 
{}

