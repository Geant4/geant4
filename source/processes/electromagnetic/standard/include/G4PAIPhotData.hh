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
// File name:     G4PAIPhotData
//
// Author:        V. Grichine based on  G4PAIModelData code for MT
//
// Creation date: 07.10.2013
//
// Modifications:
//
//
// Class Description:
//
// Implementation of PAIPhotModel internal data class.
// This class is extracted from G4PAIPhot in order to provide sharing
// of these data between threads.
//
// Internal data tables are computed for proton. 
//
// -------------------------------------------------------------------
//

#ifndef G4PAIPhotData_h
#define G4PAIPhotData_h 1

#include <vector>
#include "globals.hh"
#include "G4PAIxSection.hh"
#include "G4SandiaTable.hh"

class G4PhysicsLogVector;
class G4PhysicsTable;
class G4MaterialCutsCouple;
class G4PAIPhotModel;

class G4PAIPhotData 
{

public:

  explicit G4PAIPhotData(G4double tmin, G4double tmax, G4int verbose);

  ~G4PAIPhotData();

  void Initialise(const G4MaterialCutsCouple*, G4double cut, G4PAIPhotModel*);

  G4double DEDXPerVolume(G4int coupleIndex, G4double scaledTkin,
			 G4double cut) const;

  G4double CrossSectionPerVolume(G4int coupleIndex, G4double scaledTkin,
				 G4double tcut, G4double tmax) const;

  G4double GetPlasmonRatio( G4int coupleIndex, G4double scaledTkin ) const;

  G4double SampleAlongStepTransfer(G4int coupleIndex, G4double kinEnergy,
				   G4double scaledTkin,
				   G4double stepFactor) const;
  G4double SampleAlongStepPhotonTransfer(G4int coupleIndex, G4double kinEnergy,
				   G4double scaledTkin,
				   G4double stepFactor) const;
  G4double SampleAlongStepPlasmonTransfer(G4int coupleIndex, G4double kinEnergy,
				   G4double scaledTkin,
				   G4double stepFactor) const;

  G4double SamplePostStepTransfer(G4int coupleIndex, 
				  G4double scaledTkin) const;
  G4double SamplePostStepPhotonTransfer(G4int coupleIndex, 
				  G4double scaledTkin) const;
  G4double SamplePostStepPlasmonTransfer(G4int coupleIndex, 
				  G4double scaledTkin) const;

private:

  G4double GetEnergyTransfer(G4int coupleIndex, size_t iPlace, 
			     G4double position) const;
  G4double GetEnergyPhotonTransfer(G4int coupleIndex, size_t iPlace, 
			     G4double position) const;
  G4double GetEnergyPlasmonTransfer(G4int coupleIndex, size_t iPlace, 
			     G4double position) const;

  // hide assignment operator 
  G4PAIPhotData & operator=(const  G4PAIPhotData &right) = delete;
  G4PAIPhotData(const  G4PAIPhotData&) = delete;

  G4int                fTotBin;
  G4double             fLowestKineticEnergy;
  G4double             fHighestKineticEnergy;

  G4PhysicsLogVector*  fParticleEnergyVector;

  G4PAIxSection        fPAIxSection;
  G4SandiaTable        fSandia;

  std::vector<G4PhysicsTable*>      fPAIxscBank;
  std::vector<G4PhysicsTable*>      fPAIphotonBank;
  std::vector<G4PhysicsTable*>      fPAIplasmonBank;

  std::vector<G4PhysicsTable*>      fPAIdEdxBank;
  std::vector<G4PhysicsLogVector*>  fdEdxTable;

  std::vector<G4PhysicsLogVector*>  fdNdxCutTable;
  std::vector<G4PhysicsLogVector*>  fdNdxCutPhotonTable;
  std::vector<G4PhysicsLogVector*>  fdNdxCutPlasmonTable;

  std::vector<G4PhysicsLogVector*>  fdEdxCutTable;

};

#endif







