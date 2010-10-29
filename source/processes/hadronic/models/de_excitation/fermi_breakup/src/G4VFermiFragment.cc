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
// $Id: G4VFermiFragment.cc,v 1.6 2010-10-29 17:35:03 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Hadronic Process: Nuclear De-excitations
// by V. Lara (Nov 1998)

#include "G4VFermiFragment.hh"
#include "G4NucleiProperties.hh"

G4VFermiFragment::G4VFermiFragment(G4int anA, G4int aZ, G4int Pol, G4double ExE):
    A(anA),
    Z(aZ),
    Polarization(Pol),
    ExcitEnergy(ExE)
{
  fragmentMass = 0.0;
  if(A > 0) { fragmentMass = G4NucleiProperties::GetNuclearMass(A, Z); }
}

G4VFermiFragment::G4VFermiFragment():
    A(0),
    Z(0),
    Polarization(0),
    ExcitEnergy(0.0)
{
  fragmentMass = 0.0;
}

G4VFermiFragment::~G4VFermiFragment()
{}




