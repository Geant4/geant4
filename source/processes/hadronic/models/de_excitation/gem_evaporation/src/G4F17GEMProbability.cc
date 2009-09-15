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
// $Id: G4F17GEMProbability.cc,v 1.6 2009-09-15 12:54:17 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4F17GEMProbability.hh"

G4F17GEMProbability::G4F17GEMProbability() :
  G4GEMProbability(17,9,5.0/2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(495.33*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(286.0*picosecond);

  ExcitEnergies.push_back(3104.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(19.0*keV));

  ExcitEnergies.push_back(3857.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(4.2e-3*picosecond);

  ExcitEnergies.push_back(3857.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(4.2e-3*picosecond);

  ExcitEnergies.push_back(4696.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(225.0*keV));

  ExcitEnergies.push_back(5103.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1530.0*keV));

  ExcitEnergies.push_back(5521.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(68.0*keV));

  ExcitEnergies.push_back(5672.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(40.0*keV));

  ExcitEnergies.push_back(5682.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.6*keV));

  ExcitEnergies.push_back(5817.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(180.0*keV));

  ExcitEnergies.push_back(6036.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30.0*keV));

  ExcitEnergies.push_back(6556.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(200.0*keV));

  ExcitEnergies.push_back(6699.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(3.0*keV));

  ExcitEnergies.push_back(6774.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(4.5*keV));

  ExcitEnergies.push_back(7027.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(3.8*keV));

  ExcitEnergies.push_back(7356.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(10.0*keV));

  ExcitEnergies.push_back(7479.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(795.0*keV));

  ExcitEnergies.push_back(7546.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30.0*keV));

  ExcitEnergies.push_back(7750.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(179.0*keV));

  ExcitEnergies.push_back(8075.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

  ExcitEnergies.push_back(8200.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(700.0*keV));

  ExcitEnergies.push_back(8383.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(11.0*keV));

  ExcitEnergies.push_back(8416.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(45.0*keV));

  ExcitEnergies.push_back(8750.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(170.0*keV));

  ExcitEnergies.push_back(8760.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(90.0*keV));

  ExcitEnergies.push_back(8970.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(165.0*keV));

  ExcitEnergies.push_back(9270.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(140.0*keV));

  ExcitEnergies.push_back(9910.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(90.0*keV));

  ExcitEnergies.push_back(10040.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(280.0*keV));

  ExcitEnergies.push_back(10400.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(160.0*keV));

  ExcitEnergies.push_back(10499.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(165.0*keV));

  ExcitEnergies.push_back(10.91E3*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(560.0*keV));

  ExcitEnergies.push_back(11193.1*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.20*keV));

  ExcitEnergies.push_back(12250.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

  ExcitEnergies.push_back(12355.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(190.0*keV));

  ExcitEnergies.push_back(12500.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(660.0*keV));

  ExcitEnergies.push_back(12550.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.83*keV));

  ExcitEnergies.push_back(13061.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.0*keV));

  ExcitEnergies.push_back(13080.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.0*keV));

  ExcitEnergies.push_back(13130.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(520.0*keV));

  ExcitEnergies.push_back(13781.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(12.0*keV));

  ExcitEnergies.push_back(14000.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(260.0*keV));

  ExcitEnergies.push_back(14176.0*keV);
  ExcitSpins.push_back(3.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(30.0*keV));

  ExcitEnergies.push_back(14304.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(19.3*keV));

  ExcitEnergies.push_back(14380.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(610.0*keV));

  ExcitEnergies.push_back(14.71E3*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(470.0*keV));

  ExcitEnergies.push_back(14809.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(190.0*keV));

  ExcitEnergies.push_back(17100.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1500.0*keV));

  ExcitEnergies.push_back(19420.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(300.0*keV));

  ExcitEnergies.push_back(20250.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(350.0*keV));

  ExcitEnergies.push_back(20900.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(600.0*keV));

  ExcitEnergies.push_back(21010.0*keV);
  ExcitSpins.push_back(1.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(280.0*keV));

  ExcitEnergies.push_back(21800.0*keV);
  ExcitSpins.push_back(9.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(400.0*keV));

  ExcitEnergies.push_back(22700.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(600.0*keV));

  ExcitEnergies.push_back(23800.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(600.0*keV));

  ExcitEnergies.push_back(25400.0*keV);
  ExcitSpins.push_back(7.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1500.0*keV));

  ExcitEnergies.push_back(27200.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1500.0*keV));

  ExcitEnergies.push_back(28900.0*keV);
  ExcitSpins.push_back(5.0/2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2000.0*keV));

}


G4F17GEMProbability::G4F17GEMProbability(const G4F17GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4F17GEMProbability::copy_constructor meant to not be accessable");
}




const G4F17GEMProbability & G4F17GEMProbability::
operator=(const G4F17GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4F17GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4F17GEMProbability::operator==(const G4F17GEMProbability &) const
{
  return false;
}

G4bool G4F17GEMProbability::operator!=(const G4F17GEMProbability &) const
{
  return true;
}



