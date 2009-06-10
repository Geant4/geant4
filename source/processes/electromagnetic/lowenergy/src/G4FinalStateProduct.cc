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
// $Id: G4FinalStateProduct.cc,v 1.7 2009-06-10 13:32:36 mantero Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// Reference: TNS Geant4-DNA paper
// Reference for implementation model: NIM. 155, pp. 145-156, 1978

// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Geant4-DNA dummy final state for test purpose
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------

#include "G4FinalStateProduct.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"

G4FinalStateProduct::G4FinalStateProduct() : killStatus(false), doNotDepositStatus(false), isModified(false), localEnergyDeposit(0.), modifiedEnergy(0)
{
  // empty
}

G4FinalStateProduct::~G4FinalStateProduct()
{
  // Ownership of DynamicParticle secondaries produced is transferred
  // to ParticleChange in process 
}
 
void G4FinalStateProduct::Clear()
{
  // Reset object status
  killStatus = false;
  doNotDepositStatus = false;
  isModified = false;
  localEnergyDeposit = 0.;
  modifiedEnergy = 0.;
  secondaries.clear();
  modifiedDirection.setX(0.);
  modifiedDirection.setY(0.);
  modifiedDirection.setZ(0.);
}

void G4FinalStateProduct::AddSecondary(G4DynamicParticle* particle)
{
  secondaries.push_back(particle);
}
  
void G4FinalStateProduct::AddEnergyDeposit(G4double energy)
{
  localEnergyDeposit += energy;
}
  
G4int G4FinalStateProduct::NumberOfSecondaries() const
{
  size_t n = secondaries.size();
  return n;
}   
 
const std::vector<G4DynamicParticle*>& G4FinalStateProduct::GetSecondaries() const
{
  return secondaries;
}
  
void G4FinalStateProduct::DoNotDepositEnergy()
{
  doNotDepositStatus = true;
}

void G4FinalStateProduct::KillPrimaryParticle()
{
  
  // ---- MGP ---- To be added: Handle local energy deposit here
  killStatus = true;
}

void G4FinalStateProduct::ModifyPrimaryParticle(G4double dirX, G4double dirY, G4double dirZ, G4double energy)
{
  isModified = true;
  modifiedEnergy = energy;
  modifiedDirection.set(dirX,dirY,dirZ);
  modifiedDirection = modifiedDirection.unit();
}

void G4FinalStateProduct::ModifyPrimaryParticle(const G4ThreeVector& direction, G4double energy)
{
  isModified = true;
  modifiedEnergy = energy;
  modifiedDirection = direction.unit();
}
