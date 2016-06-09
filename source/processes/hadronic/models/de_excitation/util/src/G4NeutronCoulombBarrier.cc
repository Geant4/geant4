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
// $Id: G4NeutronCoulombBarrier.cc,v 1.3 2005/06/04 13:29:20 jwellisc Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4NeutronCoulombBarrier.hh"

G4NeutronCoulombBarrier::G4NeutronCoulombBarrier(const G4NeutronCoulombBarrier & ) : G4CoulombBarrier()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4NeutronCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4NeutronCoulombBarrier & G4NeutronCoulombBarrier::operator=(const G4NeutronCoulombBarrier & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4NeutronCoulombBarrier::operator= meant to not be accessable.");
    return *this;
}

G4bool G4NeutronCoulombBarrier::operator==(const G4NeutronCoulombBarrier & ) const 
{
    return false;
}

G4bool G4NeutronCoulombBarrier::operator!=(const G4NeutronCoulombBarrier & ) const 
{
    return true;
}

