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
// $Id: G4Ne18GEMProbability.cc,v 1.1 2002/06/06 18:02:27 larazb Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4Ne18GEMProbability.hh"

G4Ne18GEMProbability::G4Ne18GEMProbability() :
  G4GEMProbability(18,10,0.0) // A,Z,Spin
{

  ExcitEnergies.push_back(1887.3*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.34*picosecond);

  ExcitEnergies.push_back(3376.2*keV);
  ExcitSpins.push_back(4.0);
  ExcitLifetimes.push_back(3.0*picosecond);

  ExcitEnergies.push_back(3576.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(1.4*picosecond);

  ExcitEnergies.push_back(3616.4*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(0.040*picosecond);

  ExcitEnergies.push_back(4510.0*keV);
  ExcitSpins.push_back(1.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(40.0*keV));

  ExcitEnergies.push_back(4580.0*keV);
  ExcitSpins.push_back(0.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(40.0*keV));

  ExcitEnergies.push_back(7062.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(180.0*keV));

  ExcitEnergies.push_back(7915.0*keV);
  ExcitSpins.push_back(2.0);
  ExcitLifetimes.push_back(hbar_Planck*log(2.0)/(50.0*keV));

}


G4Ne18GEMProbability::G4Ne18GEMProbability(const G4Ne18GEMProbability &right)
{
  G4Exception("G4Ne18GEMProbability::copy_constructor meant to not be accessable");
}




const G4Ne18GEMProbability & G4Ne18GEMProbability::
operator=(const G4Ne18GEMProbability &right)
{
  G4Exception("G4Ne18GEMProbability::operator= meant to not be accessable");
  return *this;
}


G4bool G4Ne18GEMProbability::operator==(const G4Ne18GEMProbability &right) const
{
  return false;
}

G4bool G4Ne18GEMProbability::operator!=(const G4Ne18GEMProbability &right) const
{
  return true;
}



