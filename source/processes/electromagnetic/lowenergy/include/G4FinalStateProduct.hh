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
// $Id: G4FinalStateProduct.hh,v 1.1 2007-10-07 12:52:18 pia Exp $
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
  
  G4int NumberOfSecondaries() const;
  
  G4double GetEnergyDeposit() const;
  
  const std::vector<G4DynamicParticle*>& GetSecondaries() const;
 
  void KillIncidentParticle();
 
  G4bool GetKillParticleStatus() const;
  
  
  // protected: 
  
  // Copy constructor and assignment operator to be added here
  
private:
  
  G4bool killStatus;
  G4double localEnergyDeposit;
  std::vector<G4DynamicParticle*> secondaries;
  
};

#endif
