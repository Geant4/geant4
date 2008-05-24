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
// $Id: G4NeutronEvaporationProbability.cc,v 1.10 2008-05-24 16:34:33 ahoward Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1999)
//
//J. M. Quesada (Apr. 2008) unused items have been removed  (ExcitEnergies, ExcitSpins)


#include "G4NeutronEvaporationProbability.hh"

G4NeutronEvaporationProbability::G4NeutronEvaporationProbability() :
    G4EvaporationProbability(1,0,2,&theCoulombBarrier) // A,Z,Gamma
{
  
}


G4NeutronEvaporationProbability::G4NeutronEvaporationProbability(const G4NeutronEvaporationProbability &) : G4EvaporationProbability()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4NeutronEvaporationProbability::copy_constructor meant to not be accessable");
}




const G4NeutronEvaporationProbability & G4NeutronEvaporationProbability::
operator=(const G4NeutronEvaporationProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4NeutronEvaporationProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4NeutronEvaporationProbability::operator==(const G4NeutronEvaporationProbability &) const
{
    return false;
}

G4bool G4NeutronEvaporationProbability::operator!=(const G4NeutronEvaporationProbability &) const
{
    return true;
}
