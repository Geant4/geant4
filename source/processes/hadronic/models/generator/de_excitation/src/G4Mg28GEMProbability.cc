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
// $Id: G4Mg28GEMProbability.cc,v 1.1 2002/06/06 18:02:15 larazb Exp $
// GEANT4 tag $Name: geant4-04-01-patch-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Mg28GEMProbability.hh"

G4Mg28GEMProbability::G4Mg28GEMProbability() :
  G4GEMProbability(28,12,0.0) // A,Z,Spin
{

    ExcitEnergies.push_back(1473.4*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(1.7*picosecond);

    ExcitEnergies.push_back(3862.7*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(800.0e-3*picosecond);

    ExcitEnergies.push_back(4020.2*keV);
    ExcitSpins.push_back(4.0);
    ExcitLifetimes.push_back(150.0e-3*picosecond);

    ExcitEnergies.push_back(4557.0*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(40.0e-3*picosecond);

    ExcitEnergies.push_back(4878.6*keV);
    ExcitSpins.push_back(2.0);
    ExcitLifetimes.push_back(120.0e-3*picosecond);

    ExcitEnergies.push_back(5171.8*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(170.0e-3*picosecond);

    ExcitEnergies.push_back(5192.7*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(30.0e-3*picosecond);

    ExcitEnergies.push_back(5271.7*keV);
    ExcitSpins.push_back(1.0);
    ExcitLifetimes.push_back(150.0e-3*picosecond);

    ExcitEnergies.push_back(5702.3*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(300.0e-3*picosecond);

}


G4Mg28GEMProbability::G4Mg28GEMProbability(const G4Mg28GEMProbability &right)
{
  G4Exception("G4Mg28GEMProbability::copy_constructor meant to not be accessable");
}




const G4Mg28GEMProbability & G4Mg28GEMProbability::
operator=(const G4Mg28GEMProbability &right)
{
  G4Exception("G4Mg28GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Mg28GEMProbability::operator==(const G4Mg28GEMProbability &right) const
{
  return false;
}

G4bool G4Mg28GEMProbability::operator!=(const G4Mg28GEMProbability &right) const
{
  return true;
}



