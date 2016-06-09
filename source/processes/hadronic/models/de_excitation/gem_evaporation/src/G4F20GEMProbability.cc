//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4F20GEMProbability.cc,v 1.4 2005/06/04 13:25:25 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4F20GEMProbability.hh"

G4F20GEMProbability::G4F20GEMProbability() :
  G4GEMProbability(20,9,2.0) // A,Z,Spin
{

  ExcitEnergies.push_back(655.95*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(0.28*picosecond);

  ExcitEnergies.push_back(822.9*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(55*picosecond);

  ExcitEnergies.push_back(983.8*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(1.1*picosecond);

  ExcitEnergies.push_back(1056.93*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(31.0e-3*picosecond);

  ExcitEnergies.push_back(1309.22*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.62*picosecond);

  ExcitEnergies.push_back(1843.4*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(2043.9*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(26.0e-3*picosecond);

  ExcitEnergies.push_back(2194.6*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(8.0e-3*picosecond);

  ExcitEnergies.push_back(2966.2*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(42.0e-3*picosecond);

  ExcitEnergies.push_back(3488.4*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(30.0e-3*picosecond);

  ExcitEnergies.push_back(3525.9*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(3587.1*keV);
  ExcitSpins.push_back(3.0);
  ExcitLifetimes.push_back(21.0e-3*picosecond);

  ExcitEnergies.push_back(6627.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.29*keV));

  ExcitEnergies.push_back(6648.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.62*keV));

  ExcitEnergies.push_back(6685.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(3.80*keV));

  ExcitEnergies.push_back(6692.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(5.23*keV));

  ExcitEnergies.push_back(6696.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.05*keV));

  ExcitEnergies.push_back(6699.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(2.85*keV));

  ExcitEnergies.push_back(6709.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.14*keV));

  ExcitEnergies.push_back(6717.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.95*keV));

  ExcitEnergies.push_back(6791.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.9*keV));

  ExcitEnergies.push_back(6835.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.7*keV));

  ExcitEnergies.push_back(6837.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(0.4*keV));

  ExcitEnergies.push_back(6856.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(1.3*keV));

  ExcitEnergies.push_back(6858.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(19.0*keV));

  ExcitEnergies.push_back(7005.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(24.0*keV));

  ExcitEnergies.push_back(7076.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(24.0*keV));

  ExcitEnergies.push_back(7171.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(14.0*keV));

  ExcitEnergies.push_back(7311.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(33.0*keV));

  ExcitEnergies.push_back(7355.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(19.0*keV));

  ExcitEnergies.push_back(7410.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(10.0*keV));

  ExcitEnergies.push_back(7489.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(57.0*keV));

  ExcitEnergies.push_back(7503.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(85.0*keV));

  ExcitEnergies.push_back(7670.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(60.0*keV));

  ExcitEnergies.push_back(7800.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(100.0*keV));

  ExcitEnergies.push_back(8150.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(190.0*keV));

  ExcitEnergies.push_back(10228.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(200.0*keV));

  ExcitEnergies.push_back(10641.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(60.0*keV));

  ExcitEnergies.push_back(10807.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(330.0*keV));

}


G4F20GEMProbability::G4F20GEMProbability(const G4F20GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4F20GEMProbability::copy_constructor meant to not be accessable");
}




const G4F20GEMProbability & G4F20GEMProbability::
operator=(const G4F20GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4F20GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4F20GEMProbability::operator==(const G4F20GEMProbability &) const
{
  return false;
}

G4bool G4F20GEMProbability::operator!=(const G4F20GEMProbability &) const
{
  return true;
}



