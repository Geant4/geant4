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
// $Id: G4LossTableBuilder.hh 96088 2016-03-14 16:03:38Z gcosmo $
//
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

class G4VEmModel;
class G4ParticleDefinition;
class G4EmParameters;

class G4LossTableBuilder
{

public:

  G4LossTableBuilder();

  virtual ~G4LossTableBuilder();

  // build sum of all energy loss processes
  void BuildDEDXTable(G4PhysicsTable* dedxTable, 
		      const std::vector<G4PhysicsTable*>&);

  // build range
  void BuildRangeTable(const G4PhysicsTable* dedxTable, 
		       G4PhysicsTable* rangeTable,
		       G4bool isIonisation = false);

  // build inverse range
  void BuildInverseRangeTable(const G4PhysicsTable* rangeTable,
			      G4PhysicsTable* invRangeTable,
			      G4bool isIonisation = false);

  // build a table requested by any model class
  G4PhysicsTable* BuildTableForModel(G4PhysicsTable* table, 
				     G4VEmModel* model,
				     const G4ParticleDefinition*,
				     G4double emin, G4double emax, 
				     G4bool spline);

  // initialise base materials
  void InitialiseBaseMaterials(G4PhysicsTable* table);


  // access methods
  inline const std::vector<G4int>* GetCoupleIndexes();

  inline const std::vector<G4double>* GetDensityFactors();

  inline G4bool GetFlag(size_t idx) const;

  inline void SetSplineFlag(G4bool flag);

  inline void SetInitialisationFlag(G4bool flag);
 
private:

  void InitialiseCouples();

  G4LossTableBuilder & operator=(const  G4LossTableBuilder &right) = delete;
  G4LossTableBuilder(const  G4LossTableBuilder&) = delete;

  G4EmParameters*        theParameters;

  G4bool splineFlag;
  G4bool isInitialized;

  std::vector<G4double>* theDensityFactor;
  std::vector<G4int>*    theDensityIdx;
  std::vector<G4bool>*   theFlag;

};

inline const std::vector<G4int>* 
G4LossTableBuilder::GetCoupleIndexes()
{
  if(theDensityIdx->size() == 0) { InitialiseCouples(); }
  return theDensityIdx;
}

inline const std::vector<G4double>* 
G4LossTableBuilder::GetDensityFactors()
{
  if(theDensityIdx->size() == 0) { InitialiseCouples(); }
  return theDensityFactor;
}

inline G4bool G4LossTableBuilder::GetFlag(size_t idx) const
{
  return (*theFlag)[idx];
}

inline void G4LossTableBuilder::SetSplineFlag(G4bool flag)
{
  splineFlag = flag;
}

inline void G4LossTableBuilder::SetInitialisationFlag(G4bool flag)
{
  isInitialized = flag;
}

//....oooOO0OOooo.......oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
