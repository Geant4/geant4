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
// $Id: G4Mg27GEMProbability.cc,v 1.1 2002/06/06 18:02:14 larazb Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Mg27GEMProbability.hh"

G4Mg27GEMProbability::G4Mg27GEMProbability() :
  G4GEMProbability(27,12,1.0/2.0) // A,Z,Spin
{

    ExcitEnergies.push_back(984.66*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(1.4*picosecond);

    ExcitEnergies.push_back(1698.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(1.2*picosecond);

    ExcitEnergies.push_back(1940.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(1.1*picosecond);

    ExcitEnergies.push_back(3109.4*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(100.0e-3*picosecond);

    ExcitEnergies.push_back(3426.9*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(100.0e-3*picosecond);

    ExcitEnergies.push_back(3475.3*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(3490.7*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(15.0e-3*picosecond);

    ExcitEnergies.push_back(3559.2*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(3760.4*keV);
    ExcitSpins.push_back(7.0/2.0);
    ExcitLifetimes.push_back(610.0e-3*picosecond);

    ExcitEnergies.push_back(3785.9*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(25.0e-3*picosecond);

    ExcitEnergies.push_back(3884.0*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(700.0e-3*picosecond);

    ExcitEnergies.push_back(4149.8*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(4398.2*keV);
    ExcitSpins.push_back(5.0/2.0);
    ExcitLifetimes.push_back(65.0e-3*picosecond);

    ExcitEnergies.push_back(4552.8*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(15.0e-3*picosecond);

    ExcitEnergies.push_back(4827.3*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(4992.3*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(5028.0*keV);
    ExcitSpins.push_back(1.0/2.0);
    ExcitLifetimes.push_back(40.0e-3*picosecond);

    ExcitEnergies.push_back(5172.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(15.0e-3*picosecond);

    ExcitEnergies.push_back(5372.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(25.0e-3*picosecond);

    ExcitEnergies.push_back(5422.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(5627.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

    ExcitEnergies.push_back(5764.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(25.0e-3*picosecond);

    ExcitEnergies.push_back(5821.0*keV);
    ExcitSpins.push_back(3.0/2.0);
    ExcitLifetimes.push_back(10.0e-3*picosecond);

}


G4Mg27GEMProbability::G4Mg27GEMProbability(const G4Mg27GEMProbability &right)
{
  G4Exception("G4Mg27GEMProbability::copy_constructor meant to not be accessable");
}




const G4Mg27GEMProbability & G4Mg27GEMProbability::
operator=(const G4Mg27GEMProbability &right)
{
  G4Exception("G4Mg27GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Mg27GEMProbability::operator==(const G4Mg27GEMProbability &right) const
{
  return false;
}

G4bool G4Mg27GEMProbability::operator!=(const G4Mg27GEMProbability &right) const
{
  return true;
}



