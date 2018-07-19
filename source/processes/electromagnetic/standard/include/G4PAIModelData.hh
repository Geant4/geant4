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
// $Id: G4PAIModelData.hh 72008 2013-07-03 08:46:39Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4PAIModelData
//
// Author:        V. Ivanchenko based on V.Grichine code of G4PAIModel
//
// Creation date: 16.08.2013
//
// Modifications:
//
// 04.10.13 V. Grichine add cut of dE/dx, redirect <dE/dx> to   std::vector<G4PhysicsLogVector*>  fdEdxTable;
//
//
// Class Description:
//
// Implementation of PAI model internal data class.
// This class is extracted from G4PAIModel in order to provide sharing
// of these data between threads.
//
// Internal data tables are computed for proton. 

// -------------------------------------------------------------------
//

#ifndef G4PAIModelData_h
#define G4PAIModelData_h 1

#include <vector>
#include "globals.hh"
#include "G4PAIySection.hh"
#include "G4SandiaTable.hh"

class G4PhysicsLogVector;
class G4PhysicsTable;
class G4MaterialCutsCouple;
class G4PAIModel;

class G4PAIModelData 
{

public:

  explicit G4PAIModelData(G4double tmin, G4double tmax, G4int verbose);

  ~G4PAIModelData();

  void Initialise(const G4MaterialCutsCouple*, G4PAIModel*);

  G4double DEDXPerVolume(G4int coupleIndex, G4double scaledTkin,
			 G4double cut) const;

  G4double CrossSectionPerVolume(G4int coupleIndex, G4double scaledTkin,
				 G4double tcut, G4double tmax) const;

  G4double SampleAlongStepTransfer(G4int coupleIndex, G4double kinEnergy,
				   G4double scaledTkin,
				   G4double tmax,
				   G4double stepFactor) const;

  G4double SamplePostStepTransfer(G4int coupleIndex, 
				  G4double scaledTkin, 
				  G4double tmin, G4double tmax) const;

private:

  G4double GetEnergyTransfer(G4int coupleIndex, size_t iPlace, 
			     G4double position) const;

  // hide assignment operator 
  G4PAIModelData & operator=(const  G4PAIModelData &right) = delete;
  G4PAIModelData(const  G4PAIModelData&) = delete;

  G4int                fTotBin;
  G4double             fLowestKineticEnergy;
  G4double             fHighestKineticEnergy;

  G4PhysicsLogVector*  fParticleEnergyVector;

  G4PAIySection        fPAIySection;
  G4SandiaTable        fSandia;

  std::vector<G4PhysicsTable*>      fPAIxscBank;
  std::vector<G4PhysicsTable*>      fPAIdEdxBank;
  std::vector<G4PhysicsLogVector*>  fdEdxTable;
};

#endif

