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
/// \file radioactivedecay/rdecay02/src/exrdmEnergyDeposition.cc
/// \brief Implementation of the exrdmEnergyDeposition class
//
#include "exrdmEnergyDeposition.hh"


//
// Default constructor
//
exrdmEnergyDeposition::exrdmEnergyDeposition()
: fEnergy(0.), fTime(0.),fWeight(1.)
{;}
//
// Specific constructor
//
exrdmEnergyDeposition::exrdmEnergyDeposition( G4double energy,
                                    G4double time,
                                    G4double weight )
  : fEnergy(energy),
    fTime(time),
    fWeight(weight)
{;}


//
// Copy constructor
//
exrdmEnergyDeposition::exrdmEnergyDeposition(
                                                   const exrdmEnergyDeposition &right )
  : fEnergy(right.fEnergy),
    fTime(right.fTime),
    fWeight(right.fWeight)
{;}

//
// Destructor
//
exrdmEnergyDeposition::~exrdmEnergyDeposition() {;}

//
// Equivalence operator
//
G4bool exrdmEnergyDeposition::operator==
                                          ( const exrdmEnergyDeposition &right ) const
{
  return fTime == right.fTime;
}

//
// Order operators
//
G4bool exrdmEnergyDeposition::operator<
                                  ( const exrdmEnergyDeposition &right ) const
{
  return fTime < right.fTime;
}

G4bool exrdmEnergyDeposition::operator<=
                                  ( const exrdmEnergyDeposition &right ) const
{
  return fTime <= right.fTime;
}

