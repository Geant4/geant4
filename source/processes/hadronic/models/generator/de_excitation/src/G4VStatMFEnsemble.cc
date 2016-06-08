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
// $Id: G4VStatMFEnsemble.cc,v 1.6 2002/12/12 19:17:24 gunter Exp $
// GEANT4 tag $Name: geant4-05-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4VStatMFEnsemble.hh"


// Copy constructor
G4VStatMFEnsemble::G4VStatMFEnsemble(const G4VStatMFEnsemble & right)
{
    G4Exception("G4VStatMFEnsemble::copy_constructor meant to not be accessable");
}

// Operators

G4VStatMFEnsemble & G4VStatMFEnsemble::
operator=(const G4VStatMFEnsemble & right)
{
    G4Exception("G4VStatMFEnsemble::operator= meant to not be accessable");
    return *this;
}


G4bool G4VStatMFEnsemble::operator==(const G4VStatMFEnsemble & right) const
{
    G4Exception("G4VStatMFEnsemble::operator== meant to not be accessable");
    return false;
}
 

G4bool G4VStatMFEnsemble::operator!=(const G4VStatMFEnsemble & right) const
{
    G4Exception("G4VStatMFEnsemble::operator!= meant to not be accessable");
    return true;
}

