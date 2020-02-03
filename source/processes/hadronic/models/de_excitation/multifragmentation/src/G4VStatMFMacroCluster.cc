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
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4VStatMFMacroCluster.hh"



// Copy constructor
G4VStatMFMacroCluster::G4VStatMFMacroCluster(const G4VStatMFMacroCluster & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFMacroCluster::copy_constructor meant to not be accessible");
}

// Operators

G4VStatMFMacroCluster & G4VStatMFMacroCluster::
operator=(const G4VStatMFMacroCluster & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFMacroCluster::operator= meant to not be accessible");
    return *this;
}


G4bool G4VStatMFMacroCluster::operator==(const G4VStatMFMacroCluster & ) const
{
//	throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFMacroCluster::operator== meant to not be accessible");
    return false;
}
 

G4bool G4VStatMFMacroCluster::operator!=(const G4VStatMFMacroCluster & ) const
{
//	throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFMacroCluster::operator!= meant to not be accessible");
    return true;
}


G4double G4VStatMFMacroCluster::CalcInvLevelDensity(void)
{
    // Calculate Inverse Density Level
    // Epsilon0*(1 + 3 /(Af - 1))
    if (theA == 1) return 0.0;
    else return
	     G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(static_cast<G4double>(theA-1)));
}

