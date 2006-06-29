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
// $Id: G4StableFermiFragment.cc,v 1.5 2006-06-29 20:13:09 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4StableFermiFragment.hh"
#include "G4HadronicException.hh"

G4StableFermiFragment::G4StableFermiFragment()
{
}

G4StableFermiFragment::G4StableFermiFragment(const G4StableFermiFragment &) : G4VFermiFragment()
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StableFermiFragment::copy_constructor meant to not be accessable");
}


G4StableFermiFragment::~G4StableFermiFragment()
{
}


const G4StableFermiFragment & G4StableFermiFragment::operator=(const G4StableFermiFragment &)
{
    throw G4HadronicException(__FILE__, __LINE__, "G4StableFermiFragment::operator= meant to not be accessable");
    return *this;
}


G4bool G4StableFermiFragment::operator==(const G4StableFermiFragment &) const
{
    return false;
}

G4bool G4StableFermiFragment::operator!=(const G4StableFermiFragment &) const
{
    return true;
}



G4FragmentVector * G4StableFermiFragment::GetFragment(const G4LorentzVector & aMomentum) const
{
    G4FragmentVector * theResult = new G4FragmentVector;

    theResult->push_back(new G4Fragment(A,Z,aMomentum));
  
    return theResult;
}
