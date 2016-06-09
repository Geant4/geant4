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
// $Id: G4B8GEMProbability.cc,v 1.4 2005/06/04 13:25:25 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//


#include "G4B8GEMProbability.hh"

G4B8GEMProbability::G4B8GEMProbability() :
  G4GEMProbability(8,5,2.0) // A,Z,Spin
{
    ExcitEnergies.push_back(778.08*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(40.0*keV));
    
    ExcitEnergies.push_back(2320.0*keV);
    ExcitSpins.push_back(3.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(350.0*keV));
    
    ExcitEnergies.push_back(10619.0*keV);
    ExcitSpins.push_back(0.0);
    ExcitLifetimes.push_back(hbar_Planck*std::log(2.0)/(60.0*keV));
    
    SetExcitationEnergiesPtr(&ExcitEnergies);
    SetExcitationSpinsPtr(&ExcitSpins);
    SetExcitationLifetimesPtr(&ExcitLifetimes);
}


G4B8GEMProbability::G4B8GEMProbability(const G4B8GEMProbability &) : G4GEMProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4B8GEMProbability::copy_constructor meant to not be accessable");
}




const G4B8GEMProbability & G4B8GEMProbability::
operator=(const G4B8GEMProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4B8GEMProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4B8GEMProbability::operator==(const G4B8GEMProbability &) const
{
    return false;
}

G4bool G4B8GEMProbability::operator!=(const G4B8GEMProbability &) const
{
    return true;
}

