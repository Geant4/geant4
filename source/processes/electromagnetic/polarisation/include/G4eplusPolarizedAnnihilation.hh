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
// $Id: G4eplusPolarizedAnnihilation.hh 96114 2016-03-16 18:51:33Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eplusPolarizedAnnihilation
//
// Author:        A. Schaelicke on base of Vladimir Ivanchenko / Michel Maire code
//
// Creation date: 02.07.2006
//
// Modifications:
// 26-07-06 modified cross section  (P. Starovoitov)
// 21-08-06 interface updated   (A. Schaelicke)
// 11-06-07, add PostStepGetPhysicalInteractionLength (A.Schalicke)
//
//
// Class Description:
//
// Polarized process of e+ annihilation into 2 gammas
//

// -------------------------------------------------------------------
//

#ifndef G4eplusPolarizedAnnihilation_h
#define G4eplusPolarizedAnnihilation_h 1

#include "G4eplusAnnihilation.hh"
#include "G4Positron.hh"
#include "G4VEmModel.hh"


class G4PolarizedAnnihilationModel;

class G4eplusPolarizedAnnihilation : public G4eplusAnnihilation
{

public:

  explicit G4eplusPolarizedAnnihilation(const G4String& name = "pol-annihil");

  virtual ~G4eplusPolarizedAnnihilation();

  // Print out of the class parameters
  virtual void PrintInfo() override;

  virtual G4double GetMeanFreePath(const G4Track& track,
				   G4double previousStepSize,
				   G4ForceCondition* condition) override;

  virtual G4double PostStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double   previousStepSize,
                             G4ForceCondition* condition) override;

  virtual void BuildPhysicsTable(const G4ParticleDefinition&) override;

private:

  void CleanTables();

  void BuildAsymmetryTables(const G4ParticleDefinition& part);

  G4double ComputeAsymmetry(G4double energy,
			    const G4MaterialCutsCouple* couple,
			    const G4ParticleDefinition& particle,
			    G4double cut,
			    G4double &tasm);
  
  G4double ComputeSaturationFactor(const G4Track& aTrack);

  G4eplusPolarizedAnnihilation& operator=(const G4eplusPolarizedAnnihilation &right);
  G4eplusPolarizedAnnihilation(const G4eplusPolarizedAnnihilation& );

  G4bool isInitialised;

  // for polarization:
  G4PolarizedAnnihilationModel* emModel;
  G4ThreeVector theTargetPolarization;

  G4PhysicsTable* theAsymmetryTable;          // table for cross section assym.
  G4PhysicsTable* theTransverseAsymmetryTable; // table for transverse cross section assym.
};

#endif
