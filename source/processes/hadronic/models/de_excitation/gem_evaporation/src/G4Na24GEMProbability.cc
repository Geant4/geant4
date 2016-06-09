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
// $Id: G4Na24GEMProbability.cc,v 1.2 2003/11/03 17:53:04 hpw Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Na24GEMProbability.hh"

G4Na24GEMProbability::G4Na24GEMProbability() :
  G4GEMProbability(24,11,4.0) // A,Z,Spin
{

    ExcitEnergies.push_back(472.29*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(20.18*millisecond);

    ExcitEnergies.push_back(563.29*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(36.0*picosecond);

    ExcitEnergies.push_back(1341.4*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(62.0e-3*picosecond);

    ExcitEnergies.push_back(1344.5*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(26.0e-3*picosecond);

    ExcitEnergies.push_back(1346.5*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(1040.0e-3*picosecond);

    ExcitEnergies.push_back(1512.54*keV);
    ExcitSpins.push_back(5.0);
    ExcitLifetimes.push_back(27.0e-3*picosecond);

    ExcitEnergies.push_back(1846.1*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(180.0e-3*picosecond);

    ExcitEnergies.push_back(1885.44*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(26.0e-3*picosecond);

    ExcitEnergies.push_back(2513.4*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(2562.5*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(2903.7*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(35.0e-3*picosecond);

    ExcitEnergies.push_back(2977.8*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(3216.8*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(15.0e-3*picosecond);

    ExcitEnergies.push_back(3371.9*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(13.0e-3*picosecond);

    ExcitEnergies.push_back(3413.4*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(3589.1*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(6.0e-3*picosecond);

    ExcitEnergies.push_back(3628.5*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(3656.5*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(3681.7*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(3745.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(3935.3*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(17.0e-3*picosecond);

    ExcitEnergies.push_back(3943.39*keV);
    ExcitSpins.push_back(6.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(3977.2*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(4048.2*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(69.0e-3*picosecond);

    ExcitEnergies.push_back(4186.9*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(4196.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(4207.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(23.0e-3*picosecond);

    ExcitEnergies.push_back(4441.6*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(35.0e-3*picosecond);

    ExcitEnergies.push_back(4562.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(4621.5*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(4692.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(24.0e-3*picosecond);

    ExcitEnergies.push_back(5044.7*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(28.0e-3*picosecond);

    ExcitEnergies.push_back(5059.9*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(49.0e-3*picosecond);

    ExcitEnergies.push_back(5194.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(7.0e-3*picosecond);

    ExcitEnergies.push_back(5250.0*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(49.0e-3*picosecond);

    ExcitEnergies.push_back(5339.4*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(5398.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(14.0e-3*picosecond);

    ExcitEnergies.push_back(5480.6*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(49.0e-3*picosecond);

    ExcitEnergies.push_back(5969.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(7.0e-3*picosecond);

    ExcitEnergies.push_back(6074.2*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(35.0e-3*picosecond);

}


G4Na24GEMProbability::G4Na24GEMProbability(const G4Na24GEMProbability &) : G4GEMProbability()
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Na24GEMProbability::copy_constructor meant to not be accessable");
}




const G4Na24GEMProbability & G4Na24GEMProbability::
operator=(const G4Na24GEMProbability &)
{
  throw G4HadronicException(__FILE__, __LINE__, "G4Na24GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Na24GEMProbability::operator==(const G4Na24GEMProbability &) const
{
  return false;
}

G4bool G4Na24GEMProbability::operator!=(const G4Na24GEMProbability &) const
{
  return true;
}



