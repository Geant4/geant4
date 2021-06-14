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
//---------------------------------------------------------------------------
//
// ClassName:    G4TablesForExtrapolator
//  
// Description:  This class keep dedx, range, inverse range tables 
//               for extrapolator
//
// Author:       24.10.14 V.Ivanchenko 
//
// Modification: 
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4TablesForExtrapolator_h
#define G4TablesForExtrapolator_h 1

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4DataVector.hh"
#include <vector>

class G4ParticleDefinition;
class G4ProductionCuts;
class G4MaterialCutsCouple;
class G4LossTableBuilder;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

enum ExtTableType
{
  fDedxElectron = 0,
  fDedxPositron,
  fDedxProton,
  fDedxMuon,
  fRangeElectron,
  fRangePositron,
  fRangeProton,
  fRangeMuon,
  fInvRangeElectron,
  fInvRangePositron,
  fInvRangeProton,
  fInvRangeMuon,
  fMscElectron
};

class G4TablesForExtrapolator 
{
public:

  explicit G4TablesForExtrapolator(G4int verb, G4int bins, G4double e1, G4double e2);

  ~G4TablesForExtrapolator();

  const G4PhysicsTable* GetPhysicsTable(ExtTableType type) const; 

  void Initialisation();

  // hide assignment operator
  G4TablesForExtrapolator & operator=
  (const G4TablesForExtrapolator &right) = delete;
  G4TablesForExtrapolator(const G4TablesForExtrapolator&) = delete;

private:

  G4PhysicsTable* PrepareTable(G4PhysicsTable*);

  void ComputeElectronDEDX(const G4ParticleDefinition* part, 
			   G4PhysicsTable* table); 

  void ComputeMuonDEDX(const G4ParticleDefinition* part, 
		       G4PhysicsTable* table); 

  void ComputeProtonDEDX(const G4ParticleDefinition* part, 
			 G4PhysicsTable* table); 

  void ComputeTrasportXS(const G4ParticleDefinition* part, 
			 G4PhysicsTable* table);

  std::vector<const G4MaterialCutsCouple*> couples;
  G4DataVector cuts;

  const G4ParticleDefinition* electron;
  const G4ParticleDefinition* positron;
  const G4ParticleDefinition* muonPlus;
  const G4ParticleDefinition* muonMinus;
  const G4ParticleDefinition* proton;
  const G4ParticleDefinition* currentParticle = nullptr;

  G4LossTableBuilder* builder = nullptr;
  G4ProductionCuts* pcuts = nullptr;

  G4PhysicsTable* dedxElectron = nullptr;
  G4PhysicsTable* dedxPositron = nullptr;
  G4PhysicsTable* dedxMuon = nullptr;
  G4PhysicsTable* dedxProton = nullptr;
  G4PhysicsTable* rangeElectron = nullptr;
  G4PhysicsTable* rangePositron = nullptr;
  G4PhysicsTable* rangeMuon = nullptr;
  G4PhysicsTable* rangeProton = nullptr;
  G4PhysicsTable* invRangeElectron = nullptr;
  G4PhysicsTable* invRangePositron = nullptr;
  G4PhysicsTable* invRangeMuon = nullptr;
  G4PhysicsTable* invRangeProton = nullptr;
  G4PhysicsTable* mscElectron = nullptr;

  G4double    emin;
  G4double    emax;
  G4double    mass = 0.0;
  G4double    charge2 = 0.0;

  G4int       verbose;
  G4int       nbins;
  G4int       nmat = 0;

  G4bool      splineFlag = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

