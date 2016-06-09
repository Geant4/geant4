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
// $Id: G4VEmissionProbability.cc,v 1.6.2.1 2009/03/04 14:56:06 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-02-patch-03 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Oct 1998)
//


#include "G4VEmissionProbability.hh"
#include "G4HadronicException.hh"


G4VEmissionProbability::G4VEmissionProbability() {}
G4VEmissionProbability::~G4VEmissionProbability() {}


G4VEmissionProbability::G4VEmissionProbability(const G4VEmissionProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VEmissionProbability::copy_constructor meant to not be accessable");
}


const G4VEmissionProbability & G4VEmissionProbability::operator=(const G4VEmissionProbability &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VEmissionProbability::operator= meant to not be accessable");
    return *this;
}


G4bool G4VEmissionProbability::operator==(const G4VEmissionProbability &) const
{
    return false;
}

G4bool G4VEmissionProbability::operator!=(const G4VEmissionProbability &) const
{
    return true;
}





