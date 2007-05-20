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
// $Id: G4InuclEvaporation.cc,v 1.2 2007-05-20 20:03:21 miheikki Exp $
//
#include "G4InuclEvaporation.hh"
#include "G4HadronicException.hh"
#include <numeric>

G4InuclEvaporation::G4InuclEvaporation() {
}

G4InuclEvaporation::G4InuclEvaporation(const G4InuclEvaporation &) : G4VEvaporation() {
    throw G4HadronicException(__FILE__, __LINE__, "G4InuclEvaporation::copy_constructor meant to not be accessable.");
}


G4InuclEvaporation::~G4InuclEvaporation() {
}

const G4InuclEvaporation & G4InuclEvaporation::operator=(const G4InuclEvaporation &) {
    throw G4HadronicException(__FILE__, __LINE__, "G4InuclEvaporation::operator= meant to not be accessable.");
    return *this;
}


G4bool G4InuclEvaporation::operator==(const G4InuclEvaporation &) const {
    return false;
}

G4bool G4InuclEvaporation::operator!=(const G4InuclEvaporation &) const {
    return true;
}

G4FragmentVector * G4InuclEvaporation::BreakItUp(const G4Fragment &theNucleus) {
    G4FragmentVector * theResult = new G4FragmentVector;

    // CHECK that Excitation Energy != 0
    if (theNucleus.GetExcitationEnergy() <= 0.0) {
	theResult->push_back(new G4Fragment(theNucleus));
	return theResult;
    }

    // The residual nucleus (after evaporation of each fragment)
    G4Fragment theResidualNucleus = theNucleus;

	

    // Starts loop over evaporated particles
    for (;;) {

#ifdef DEBUG
		G4cout <<           "-----------------------------------------------------------\n"; 
#endif  


    return theResult;
    }
}




