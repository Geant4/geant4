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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4LossTableBuilder
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
// 
// Creation date: 03.01.2002
//
// Modifications: 
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 17-07-08 Added splineFlag (V.Ivanchenko)
//
// Class Description: 
//
// Provide building of dE/dx, range, and inverse range tables.

// -------------------------------------------------------------------
//

#ifndef G4LossTableBuilder_h
#define G4LossTableBuilder_h 1

#include <vector>
#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4Threading.hh"

class G4VEmModel;
class G4ParticleDefinition;
class G4EmParameters;

class G4LossTableBuilder
{

public:

  G4LossTableBuilder(G4bool master=true);

  ~G4LossTableBuilder();

  // build sum of all energy loss processes
  void BuildDEDXTable(G4PhysicsTable* dedxTable, 
		      const std::vector<G4PhysicsTable*>&);

  // build range
  void BuildRangeTable(const G4PhysicsTable* dedxTable, 
		       G4PhysicsTable* rangeTable);

  // build inverse range
  void BuildInverseRangeTable(const G4PhysicsTable* rangeTable,
			      G4PhysicsTable* invRangeTable);

  // build a table requested by any model class
  G4PhysicsTable* BuildTableForModel(G4PhysicsTable* table, 
				     G4VEmModel* model,
				     const G4ParticleDefinition*,
				     G4double emin, G4double emax, 
				     G4bool spline);

  // initialise base materials
  void InitialiseBaseMaterials(const G4PhysicsTable* table=nullptr);

  // access methods
  const std::vector<G4int>* GetCoupleIndexes() const;

  const std::vector<G4double>* GetDensityFactors() const;

  G4bool GetFlag(size_t idx);

  G4bool GetBaseMaterialFlag();

  inline void SetSplineFlag(G4bool flag);

  inline void SetInitialisationFlag(G4bool flag);

  inline void SetBaseMaterialActive(G4bool flag);

  G4LossTableBuilder & operator=(const  G4LossTableBuilder &right) = delete;
  G4LossTableBuilder(const  G4LossTableBuilder&) = delete;
 
private:

  G4EmParameters* theParameters;

  G4bool splineFlag = true;
  G4bool isInitialized = false;
  G4bool baseMatFlag = false;
  G4bool isBaseMatActive = true;
  G4bool isMaster;

  static std::vector<G4double>* theDensityFactor;
  static std::vector<G4int>*    theDensityIdx;
  static std::vector<G4bool>*   theFlag;
#ifdef G4MULTITHREADED
  static G4Mutex ltbMutex;
#endif
};

inline void G4LossTableBuilder::SetSplineFlag(G4bool flag)
{
  splineFlag = flag;
}

inline void G4LossTableBuilder::SetInitialisationFlag(G4bool flag)
{
  isInitialized = flag;
}

inline void G4LossTableBuilder::SetBaseMaterialActive(G4bool flag)
{
  isBaseMatActive = flag;
  if(!flag) {
    baseMatFlag = false;
    isInitialized = false;
  } 
}

//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
