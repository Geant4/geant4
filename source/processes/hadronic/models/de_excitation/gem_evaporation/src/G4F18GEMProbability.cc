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
// $Id: G4F18GEMProbability.cc 87017 2014-11-21 16:26:26Z gcosmo $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//

#include "G4F18GEMProbability.hh"
#include "G4SystemOfUnits.hh"

G4F18GEMProbability::G4F18GEMProbability() :
  G4GEMProbability(18,9,1.0) // A,Z,Spin
{
  ExcitEnergies.push_back(937.2*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(47.0*picosecond);

  ExcitEnergies.push_back(1041.55*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(3.0e-3*picosecond);

  ExcitEnergies.push_back(1080.54*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(21.0*picosecond);

  ExcitEnergies.push_back(1121.36*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(151.0*nanosecond);

  ExcitEnergies.push_back(1700.81*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(0.60*picosecond);

  ExcitEnergies.push_back(2100.61*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(3.0*picosecond);

  ExcitEnergies.push_back(2523.35*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.45*picosecond);

  ExcitEnergies.push_back(3061.84*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.8e-3*picosecond);

  ExcitEnergies.push_back(3133.87*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(0.26*picosecond);

  ExcitEnergies.push_back(3358.2*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(0.33*picosecond);

  ExcitEnergies.push_back(3724.19*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(0.06*picosecond);

  ExcitEnergies.push_back(3839.17*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.050*picosecond);

  ExcitEnergies.push_back(4360.15*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(0.42*picosecond);

  ExcitEnergies.push_back(5297.6*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(21.0*picosecond);

  ExcitEnergies.push_back(5603.38*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(1.2*keV));

  ExcitEnergies.push_back(5672.57*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(0.8*keV));

  ExcitEnergies.push_back(6096.4*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(2.0*keV));

  ExcitEnergies.push_back(6163.2*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(15.0*keV));

  ExcitEnergies.push_back(6242.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(0.8*keV));

  ExcitEnergies.push_back(6262.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(3.0*keV));

  ExcitEnergies.push_back(6283.2*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(7.5*keV));

  ExcitEnergies.push_back(6310.5*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(3.1*keV));

  ExcitEnergies.push_back(6385.5*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(4.5*keV));

  ExcitEnergies.push_back(6484.9*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(1.2*keV));

  ExcitEnergies.push_back(6567.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(0.8*keV));

  ExcitEnergies.push_back(6643.7*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(89.0*keV));

  ExcitEnergies.push_back(6647.4*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(1.6*keV));

  ExcitEnergies.push_back(6777.4*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(10.0*keV));

  ExcitEnergies.push_back(6803.1*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(90.0*keV));

  ExcitEnergies.push_back(7194.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(4.0*keV));

  ExcitEnergies.push_back(7206.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(45.0*keV));

  ExcitEnergies.push_back(7313.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(53.0*keV));

  ExcitEnergies.push_back(7395.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(35.0*keV));

  ExcitEnergies.push_back(7872.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(30.0*keV));

  ExcitEnergies.push_back(7950.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(9145.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(108.0*keV));

  ExcitEnergies.push_back(9260.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(30.0*keV));

  ExcitEnergies.push_back(9320.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(9700.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(370.0*keV));

  ExcitEnergies.push_back(9820.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(10210.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(90.0*keV));

  ExcitEnergies.push_back(10268.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(180.0*keV));

  ExcitEnergies.push_back(10310.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(570.0*keV));

  ExcitEnergies.push_back(10352.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(270.0*keV));

  ExcitEnergies.push_back(10420.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

  ExcitEnergies.push_back(10541.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

  ExcitEnergies.push_back(10590.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(160.0*keV));

  ExcitEnergies.push_back(10610.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(10630.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(10640.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(240.0*keV));

  ExcitEnergies.push_back(10680.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(10690.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(75.0*keV));

  ExcitEnergies.push_back(10710.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(375.0*keV));

  ExcitEnergies.push_back(10800.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(85.0*keV));

  ExcitEnergies.push_back(10910.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(240.0*keV));

  ExcitEnergies.push_back(10950.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(400.0*keV));

  ExcitEnergies.push_back(10980.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(120.0*keV));

  ExcitEnergies.push_back(10990.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(220.0*keV));

  ExcitEnergies.push_back(11030.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(11073.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

  ExcitEnergies.push_back(11110.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(1000.0*keV));

  ExcitEnergies.push_back(11130.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(35.0*keV));

  ExcitEnergies.push_back(11260.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(240.0*keV));

  ExcitEnergies.push_back(11300.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(475.0*keV));

  ExcitEnergies.push_back(11310.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(11320.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(95.0*keV));

  ExcitEnergies.push_back(11320.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(65.0*keV));

  ExcitEnergies.push_back(11510.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(125.0*keV));

  ExcitEnergies.push_back(11560.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(11660.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(240.0*keV));

  ExcitEnergies.push_back(11830.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(180.0*keV));

  ExcitEnergies.push_back(11960.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(30.0*keV));

  ExcitEnergies.push_back(11970.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(125.0*keV));

  ExcitEnergies.push_back(12120.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(35.0*keV));

  ExcitEnergies.push_back(12140.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(220.0*keV));

  ExcitEnergies.push_back(12250.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(170.0*keV));

  ExcitEnergies.push_back(12310.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(12340.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(190.0*keV));

  ExcitEnergies.push_back(12360.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(12380.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(120.0*keV));

  ExcitEnergies.push_back(12470.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(35.0*keV));

  ExcitEnergies.push_back(12480.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(120.0*keV));

  ExcitEnergies.push_back(12484.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(12526.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(12565.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(300.0*keV));

  ExcitEnergies.push_back(12630.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(210.0*keV));

  ExcitEnergies.push_back(12670.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(290.0*keV));

  ExcitEnergies.push_back(12683.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(45.0*keV));

  ExcitEnergies.push_back(12702.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(12810.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(440.0*keV));

  ExcitEnergies.push_back(12900.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(120.0*keV));

  ExcitEnergies.push_back(12990.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(650.0*keV));

  ExcitEnergies.push_back(13030.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(250.0*keV));

  ExcitEnergies.push_back(13030.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(13129.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(35.0*keV));

  ExcitEnergies.push_back(13170.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(185.0*keV));

  ExcitEnergies.push_back(13198.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(140.0*keV));

  ExcitEnergies.push_back(13200.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(255.0*keV));

  ExcitEnergies.push_back(13260.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(360.0*keV));

  ExcitEnergies.push_back(13290.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(13332.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

  ExcitEnergies.push_back(13398.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(85.0*keV));

  ExcitEnergies.push_back(13430.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(235.0*keV));

  ExcitEnergies.push_back(13440.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(13464.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(13500.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(65.0*keV));

  ExcitEnergies.push_back(13597.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(100.0*keV));

  ExcitEnergies.push_back(13600.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(40.0*keV));

  ExcitEnergies.push_back(13610.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(90.0*keV));

  ExcitEnergies.push_back(13650.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(300.0*keV));

  ExcitEnergies.push_back(13652.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(13671.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(60.0*keV));

  ExcitEnergies.push_back(13679.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(13700.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(80.0*keV));

  ExcitEnergies.push_back(13770.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(fPlanck/(205.0*keV));

  ExcitEnergies.push_back(13780.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(130.0*keV));

  ExcitEnergies.push_back(13810.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(265.0*keV));

  ExcitEnergies.push_back(13960.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(fPlanck/(70.0*keV));

  ExcitEnergies.push_back(13980.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(170.0*keV));

  ExcitEnergies.push_back(14100.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(fPlanck/(170.0*keV));

  ExcitEnergies.push_back(14162.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(14323.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(fPlanck/(95.0*keV));

  ExcitEnergies.push_back(16513.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

  ExcitEnergies.push_back(16606.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(fPlanck/(50.0*keV));

  ExcitEnergies.push_back(16723.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

  ExcitEnergies.push_back(16870.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(fPlanck/(200.0*keV));

  ExcitEnergies.push_back(16948.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(fPlanck/(150.0*keV));

}

G4F18GEMProbability::~G4F18GEMProbability() 
{}
