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
// $Id: G4VCoulombBarrier.cc,v 1.7 2009-03-04 11:05:02 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Dec 1999)


#include "G4VCoulombBarrier.hh"
#include "G4HadronicException.hh"
#include <sstream>

G4VCoulombBarrier::G4VCoulombBarrier()
  : theA(1),theZ(0)
{
}


G4VCoulombBarrier::G4VCoulombBarrier(const G4int anA, const G4int aZ)
{
    if (anA >= aZ && anA > 0) {
	theA = anA;
	theZ = aZ;
    } else {
        std::ostringstream errOs;
	errOs << "G4VCoulombBarrier::G4VCoulombBarrier: ";
	errOs << "Wrong values for ";
	errOs << "A = " << anA << " ";
	errOs << "and Z = " << aZ << G4endl;
	throw G4HadronicException(__FILE__, __LINE__, errOs.str());
    }
}


G4VCoulombBarrier::~G4VCoulombBarrier()
{
}


G4VCoulombBarrier::G4VCoulombBarrier(const G4VCoulombBarrier & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VCoulombBarrier::copy_constructor meant to not be accessable.");
}


const G4VCoulombBarrier & G4VCoulombBarrier::operator=(const G4VCoulombBarrier & )
{
    throw G4HadronicException(__FILE__, __LINE__, "G4VCoulombBarrier::operator= meant to not be accessable.");
    return *this;
}

G4bool G4VCoulombBarrier::operator==(const G4VCoulombBarrier & ) const 
{
    return false;
}

G4bool G4VCoulombBarrier::operator!=(const G4VCoulombBarrier & ) const 
{
    return true;
}

