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
// $Id: G4VStatMFMacroCluster.cc,v 1.4.2.1 2001/06/28 19:13:23 gunter Exp $
// GEANT4 tag $Name:  $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4VStatMFMacroCluster.hh"



// Copy constructor
G4VStatMFMacroCluster::G4VStatMFMacroCluster(const G4VStatMFMacroCluster & right)
{
    G4Exception("G4VStatMFMacroCluster::copy_constructor meant to not be accessable");
}

// Operators

G4VStatMFMacroCluster & G4VStatMFMacroCluster::
operator=(const G4VStatMFMacroCluster & right)
{
    G4Exception("G4VStatMFMacroCluster::operator= meant to not be accessable");
    return *this;
}


G4bool G4VStatMFMacroCluster::operator==(const G4VStatMFMacroCluster & right) const
{
//	G4Exception("G4VStatMFMacroCluster::operator== meant to not be accessable");
    return false;
}
 

G4bool G4VStatMFMacroCluster::operator!=(const G4VStatMFMacroCluster & right) const
{
//	G4Exception("G4VStatMFMacroCluster::operator!= meant to not be accessable");
    return true;
}


G4double G4VStatMFMacroCluster::CalcInvLevelDensity(void)
{
    // Calculate Inverse Density Level
    // Epsilon0*(1 + 3 /(Af - 1))
    if (theA == 1) return 0.0;
    else return
	     G4StatMFParameters::GetEpsilon0()*(1.0+3.0/(G4double(theA) - 1.0));
}

