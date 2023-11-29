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
#include "G4SDParticleFilter.hh"
#include "G4Step.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

////////////////////////////////////////////////////////////////////////////////
// class description:
//
//  This is the class of a filter to be associated with a
// sensitive detector.
//  This class filters steps by partilce definition.
//
// Created: 2005-11-14  Tsukasa ASO.
//
///////////////////////////////////////////////////////////////////////////////

G4SDParticleFilter::G4SDParticleFilter(G4String name)
  : G4VSDFilter(name)
{}

G4SDParticleFilter::G4SDParticleFilter(G4String name,
                                       const G4String& particleName)
  : G4VSDFilter(name)
{
  G4ParticleDefinition* pd =
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(pd == nullptr)
  {
    G4String msg = "Particle <";
    msg += particleName;
    msg += "> not found.";
    G4Exception("G4SDParticleFilter::G4SDParticleFilter", "DetPS0101",
                FatalException, msg);
  }
  thePdef.push_back(pd);
}

G4SDParticleFilter::G4SDParticleFilter(
  G4String name, const std::vector<G4String>& particleNames)
  : G4VSDFilter(name)
{
  for(const auto & particleName : particleNames)
  {
    G4ParticleDefinition* pd =
      G4ParticleTable::GetParticleTable()->FindParticle(particleName);
    if(pd == nullptr)
    {
      G4String msg = "Particle <";
      msg += particleName;
      msg += "> not found.";
      G4Exception("G4SDParticleFilter::G4SDParticleFilter", "DetPS0102",
                  FatalException, msg);
    }
    thePdef.push_back(pd);
  }
}

G4SDParticleFilter::G4SDParticleFilter(
  G4String name, const std::vector<G4ParticleDefinition*>& particleDef)
  : G4VSDFilter(name)
  , thePdef(particleDef)
{
  for(auto i : particleDef)
  {
    if(i == nullptr)
      G4Exception("G4SDParticleFilter::G4SDParticleFilter", "DetPS0103",
                  FatalException,
                  "NULL pointer is found in the given particleDef vector.");
  }
}

G4bool G4SDParticleFilter::Accept(const G4Step* aStep) const
{
  for(auto i : thePdef)
  {
    if(i == aStep->GetTrack()->GetDefinition())
      return true;
  }

  // Ions by Z,A
  for(size_t i = 0; i < theIonZ.size(); i++)
  {
    if(theIonZ[i] == aStep->GetTrack()->GetDefinition()->GetAtomicNumber() &&
       theIonA[i] == aStep->GetTrack()->GetDefinition()->GetAtomicMass())
    {
      return true;
    }
  }

  return false;
}

void G4SDParticleFilter::add(const G4String& particleName)
{
  G4ParticleDefinition* pd =
    G4ParticleTable::GetParticleTable()->FindParticle(particleName);
  if(pd == nullptr)
  {
    G4String msg = "Particle <";
    msg += particleName;
    msg += "> not found.";
    G4Exception("G4SDParticleFilter::add()", "DetPS0104", FatalException, msg);
  }
  for(auto & i : thePdef)
  {
    if(i == pd)
      return;
  }
  thePdef.push_back(pd);
}

void G4SDParticleFilter::addIon(G4int Z, G4int A)
{
  for(size_t i = 0; i < theIonZ.size(); i++)
  {
    if(theIonZ[i] == Z && theIonA[i] == A)
    {
      G4cout << "G4SDParticleFilter:: Ion has been already registered."
             << G4endl;
      return;
    }
  }
  theIonZ.push_back(Z);
  theIonA.push_back(A);
}

void G4SDParticleFilter::show()
{
  G4cout << "----G4SDParticleFileter particle list------" << G4endl;
  for(auto & i : thePdef)
  {
    G4cout << i->GetParticleName() << G4endl;
  }
  for(size_t i = 0; i < theIonZ.size(); i++)
  {
    G4cout << " Ion PrtclDef (" << theIonZ[i] << "," << theIonA[i] << ")"
           << G4endl;
  }
  G4cout << "-------------------------------------------" << G4endl;
}
