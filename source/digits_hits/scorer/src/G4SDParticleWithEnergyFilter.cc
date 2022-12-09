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
// G4VSensitiveDetector
#include "G4SDParticleWithEnergyFilter.hh"
#include "G4Step.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector.
//  This class filters steps by partilce definition and kinetic energy.
//
// Created: 2005-11-14  Tsukasa ASO.
//
///////////////////////////////////////////////////////////////////////////////

G4SDParticleWithEnergyFilter::G4SDParticleWithEnergyFilter(G4String name,
                                                           G4double elow,
                                                           G4double ehigh)
  : G4VSDFilter(name)
{
  fParticleFilter = new G4SDParticleFilter(name);
  fKineticFilter  = new G4SDKineticEnergyFilter(name, elow, ehigh);
}

G4SDParticleWithEnergyFilter::~G4SDParticleWithEnergyFilter()
{
  delete fParticleFilter;
  delete fKineticFilter;
}

G4bool G4SDParticleWithEnergyFilter::Accept(const G4Step* aStep) const
{
  if(!fParticleFilter->Accept(aStep))
    return false;
  if(!fKineticFilter->Accept(aStep))
    return false;
  return true;
}

void G4SDParticleWithEnergyFilter::add(const G4String& particleName)
{
  fParticleFilter->add(particleName);
}

void G4SDParticleWithEnergyFilter::SetKineticEnergy(G4double elow,
                                                    G4double ehigh)
{
  fKineticFilter->SetKineticEnergy(elow, ehigh);
}

void G4SDParticleWithEnergyFilter::show()
{
  fParticleFilter->show();
  fKineticFilter->show();
}

G4SDParticleWithEnergyFilter::G4SDParticleWithEnergyFilter(
  const G4SDParticleWithEnergyFilter& rhs)
  : G4VSDFilter(rhs.filterName)
{
  fParticleFilter = new G4SDParticleFilter(*rhs.fParticleFilter);
  fKineticFilter  = new G4SDKineticEnergyFilter(*rhs.fKineticFilter);
}

G4SDParticleWithEnergyFilter& G4SDParticleWithEnergyFilter::operator=(
  const G4SDParticleWithEnergyFilter& rhs)
{
  if(this == &rhs)
    return *this;
  G4VSDFilter::operator=(rhs);
  delete fParticleFilter;
  fParticleFilter = new G4SDParticleFilter(*(rhs.fParticleFilter));
  delete fKineticFilter;
  fKineticFilter = new G4SDKineticEnergyFilter(*(rhs.fKineticFilter));
  return *this;
}
