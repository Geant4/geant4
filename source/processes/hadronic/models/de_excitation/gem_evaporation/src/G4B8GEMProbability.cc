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
// $Id: G4B8GEMProbability.cc,v 1.6 2009-09-15 12:54:16 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

