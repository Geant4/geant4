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
// $Id: G4FinalStateProduct.hh,v 1.5 2009-01-20 07:50:28 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Contact Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// Date         Name              Modification
// 28 Apr 2007  M.G. Pia          Created in compliance with design described in TNS paper
//
// -------------------------------------------------------------------

// Class description:
// Final state product
// Reference: TNS Geant4-DNA paper
// Further documentation available from http://www.ge.infn.it/geant4/dna

// -------------------------------------------------------------------


#ifndef G4FINALSTATEPRODUCT_HH
#define G4FINALSTATEPRODUCT_HH 1
 
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include <vector>

class G4DynamicParticle;

class G4FinalStateProduct
{
public:
  
  G4FinalStateProduct();
  
  ~G4FinalStateProduct();
  
  void Clear();
  
  void AddSecondary(G4DynamicParticle* particle);
  
  void AddEnergyDeposit(G4double energy);
  
  void ModifyPrimaryParticle(G4double dirX, G4double dirY, G4double dirZ, G4double energy);
  void ModifyPrimaryParticle(const G4ThreeVector& direction, G4double energy);
  
  void DoNotDepositEnergy();
  void KillPrimaryParticle();

  G4bool PrimaryParticleIsKilled() const { return killStatus; }
  
  G4bool PrimaryParticleIsKilledAndDoNotDepositEnergy() const { return doNotDepositStatus; }
  
  G4bool PrimaryParticleIsModified() const { return isModified; }
 
  G4int NumberOfSecondaries() const;
  
  G4double GetEnergyDeposit() const { return localEnergyDeposit; }

  G4double GetModifiedEnergy() const { return modifiedEnergy; }

  const G4ThreeVector& GetModifiedDirection() const { return modifiedDirection; }
  
  const std::vector<G4DynamicParticle*>& GetSecondaries() const;
 
  
  // protected: 
  
  // Copy constructor and assignment operator to be added here
  
private:
  
  G4bool killStatus;
  G4bool doNotDepositStatus;
  G4bool isModified;
  G4double localEnergyDeposit;
  G4double modifiedEnergy;
  std::vector<G4DynamicParticle*> secondaries;
  G4ThreeVector modifiedDirection;
};

#endif
