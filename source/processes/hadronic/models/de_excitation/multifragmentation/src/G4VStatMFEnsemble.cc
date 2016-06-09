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
// $Id: G4VStatMFEnsemble.cc,v 1.4 2005/06/04 13:27:49 jwellisc Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara

#include "G4VStatMFEnsemble.hh"
#include "G4HadronicException.hh"

// Copy constructor
G4VStatMFEnsemble::G4VStatMFEnsemble(const G4VStatMFEnsemble & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFEnsemble::copy_constructor meant to not be accessable");
}

// Operators

G4VStatMFEnsemble & G4VStatMFEnsemble::
operator=(const G4VStatMFEnsemble & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFEnsemble::operator= meant to not be accessable");
    return *this;
}


G4bool G4VStatMFEnsemble::operator==(const G4VStatMFEnsemble & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFEnsemble::operator== meant to not be accessable");
    return false;
}
 

G4bool G4VStatMFEnsemble::operator!=(const G4VStatMFEnsemble & ) const
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VStatMFEnsemble::operator!= meant to not be accessable");
    return true;
}

