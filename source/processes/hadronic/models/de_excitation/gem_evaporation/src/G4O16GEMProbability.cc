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
//
// $Id: G4O16GEMProbability.cc,v 1.6 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4O16GEMProbability.hh"

G4O16GEMProbability::G4O16GEMProbability() :
  G4GEMProbability(16,8,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(6049.4*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(67.0*picosecond);

  ExcitEnergies.push_back(6130.43*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(18.4*picosecond);

  ExcitEnergies.push_back(6917.1*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(4.16E-3*picosecond);

  ExcitEnergies.push_back(7116.85*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(8.30E-3*picosecond);

  ExcitEnergies.push_back(8871.9*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(125.0E-3*picosecond);

  ExcitEnergies.push_back(9632.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(510.0*keV));

  ExcitEnergies.push_back(9847.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.9*keV));

  ExcitEnergies.push_back(10353.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(27.0*keV));

  ExcitEnergies.push_back(10952.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(6.E-3*picosecond);

  ExcitEnergies.push_back(11080.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(12.0*keV));

  ExcitEnergies.push_back(11095.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.28*keV));

  ExcitEnergies.push_back(11260.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2500.0*keV));

  ExcitEnergies.push_back(11521.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(74.0*keV));

  ExcitEnergies.push_back(11600.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(800.0*keV));

  ExcitEnergies.push_back(12053.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.5*keV));

  ExcitEnergies.push_back(12442.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(98.0*keV));

  ExcitEnergies.push_back(12530.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.8*keV));

  ExcitEnergies.push_back(12797.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(38.0*keV));

  ExcitEnergies.push_back(12968.6*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.9*keV));

  ExcitEnergies.push_back(13020.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(150.0*keV));

  ExcitEnergies.push_back(13094.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(130.0*keV));

  ExcitEnergies.push_back(13129.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(110.0*keV));

  ExcitEnergies.push_back(13254.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(21.0*keV));

  ExcitEnergies.push_back(13664.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(64.0*keV));

  ExcitEnergies.push_back(13875.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(75.0*keV));

  ExcitEnergies.push_back(13979.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(22.0*keV));

  ExcitEnergies.push_back(14000.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(4800.0*keV));

  ExcitEnergies.push_back(14032.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(185.0*keV));

  ExcitEnergies.push_back(14.1E3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(750.0*keV));

  ExcitEnergies.push_back(14400.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30.0*keV));

  ExcitEnergies.push_back(14630.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(500.0*keV));

  ExcitEnergies.push_back(14670.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(560.0*keV));

  ExcitEnergies.push_back(14815.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(67.0*keV));

  ExcitEnergies.push_back(14922.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(51.0*keV));

  ExcitEnergies.push_back(15170.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(190.0*keV));

  ExcitEnergies.push_back(15220.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(70.0*keV));

  ExcitEnergies.push_back(15250.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(650.0*keV));

  ExcitEnergies.push_back(15450.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(95.0*keV));

  ExcitEnergies.push_back(15800.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(400.0*keV));

  ExcitEnergies.push_back(15900.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(600.0*keV));

  ExcitEnergies.push_back(16214.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(96.0*keV));

  ExcitEnergies.push_back(16220.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(18.0*keV));

  ExcitEnergies.push_back(16290.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(370.0*keV));

  ExcitEnergies.push_back(16420.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(35.0*keV));

  ExcitEnergies.push_back(16.8E3*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

  ExcitEnergies.push_back(16900.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(700.0*keV));

  ExcitEnergies.push_back(16940.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(280.0*keV));

  ExcitEnergies.push_back(17000.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1400.0*keV));

  ExcitEnergies.push_back(17140.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(36.0*keV));

  ExcitEnergies.push_back(17200.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(160.0*keV));

  ExcitEnergies.push_back(17290.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(90.0*keV));

  ExcitEnergies.push_back(17550.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(165.0*keV));

  ExcitEnergies.push_back(17640.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(59.0*keV));

  ExcitEnergies.push_back(17850.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

  ExcitEnergies.push_back(17850.0*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(14.0*keV));

  ExcitEnergies.push_back(18400.0*keV);
  ExcitSpins.push_back(5.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(510.0*keV));

  ExcitEnergies.push_back(18480.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(60.0*keV));

  ExcitEnergies.push_back(18690.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(280.0*keV));

  ExcitEnergies.push_back(18940.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(16.0*keV));

  ExcitEnergies.push_back(18990.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(240.0*keV));

  ExcitEnergies.push_back(19090.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(120.0*keV));

  ExcitEnergies.push_back(19240.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(90.0*keV));

  ExcitEnergies.push_back(19340.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(50.0*keV));

  ExcitEnergies.push_back(19480.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(250.0*keV));

  ExcitEnergies.push_back(19890.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

  ExcitEnergies.push_back(20.15E3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(350.0*keV));

  ExcitEnergies.push_back(20360.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(500.0*keV));

  ExcitEnergies.push_back(20880.0*keV);
  ExcitSpins.push_back(7.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(650.0*keV));

  ExcitEnergies.push_back(20.9E3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(350.0*keV));

  ExcitEnergies.push_back(20945.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

  ExcitEnergies.push_back(21030.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(255.0*keV));

  ExcitEnergies.push_back(21820.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(400.0*keV));

  ExcitEnergies.push_back(21840.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(55.0*keV));

  ExcitEnergies.push_back(22146.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(675.0*keV));

  ExcitEnergies.push_back(22721.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(12.5*keV));

  ExcitEnergies.push_back(22870.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

  ExcitEnergies.push_back(23100.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(500.0*keV));

  ExcitEnergies.push_back(23220.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(400.0*keV));

  ExcitEnergies.push_back(23879.0*keV);
  ExcitSpins.push_back(6.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(26.0*keV));

  ExcitEnergies.push_back(24065.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(550.0*keV));

  ExcitEnergies.push_back(24522.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(50.0*keV));

  ExcitEnergies.push_back(25120.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2900.0*keV));

  ExcitEnergies.push_back(25.50E3*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1330.0*keV));

  ExcitEnergies.push_back(26300.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1200.0*keV));

  ExcitEnergies.push_back(35000.0*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(5000.0*keV));

}


G4O16GEMProbability::G4O16GEMProbability(const G4O16GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O16GEMProbability::copy_constructor meant to not be accessable");
}




const G4O16GEMProbability & G4O16GEMProbability::
operator=(const G4O16GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4O16GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4O16GEMProbability::operator==(const G4O16GEMProbability &) const
{
  return false;
}

G4bool G4O16GEMProbability::operator!=(const G4O16GEMProbability &) const
{
  return true;
}



