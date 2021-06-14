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
//---------------------------------------------------------------------------
//
// ClassName:    G4EnergyLossForExtrapolator
//  
// Description:  This class provide calculation of energy loss, fluctuation, 
//               and msc angle
//
// Author:       09.12.04 V.Ivanchenko 
//
// Modification: 
// 08-04-05 Rename Propogator -> Extrapolator
// 16-03-06 Add muon tables
// 21-03-06 Add verbosity defined in the constructor and Initialisation
//          start only when first public method is called (V.Ivanchenko)
// 03-05-06 Remove unused pointer G4Material* from number of methods (VI)
// 28-07-07 Add maxEnergyTransfer for computation of energy loss (VI)
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef G4EnergyLossForExtrapolator_h
#define G4EnergyLossForExtrapolator_h 1

#include <vector>
#include <CLHEP/Units/PhysicalConstants.h>

#include "globals.hh"
#include "G4PhysicsTable.hh"
#include "G4TablesForExtrapolator.hh"
#include "G4Log.hh"
#include "G4Threading.hh"

class G4ParticleDefinition;
class G4Material;
class G4MaterialCutsCouple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EnergyLossForExtrapolator 
{
public:

  explicit G4EnergyLossForExtrapolator(G4int verb = 1);

  ~G4EnergyLossForExtrapolator();

  void Initialisation();

  G4double ComputeDEDX(G4double kinEnergy, const G4ParticleDefinition*,
                       const G4Material*);

  G4double ComputeRange(G4double kinEnergy, const G4ParticleDefinition*,
                       const G4Material*);

  G4double ComputeEnergy(G4double range, const G4ParticleDefinition*,
                       const G4Material*);

  G4double EnergyAfterStep(G4double kinEnergy, G4double step, 
			   const G4Material*, const G4ParticleDefinition*);

  G4double EnergyBeforeStep(G4double kinEnergy, G4double step, 
	        	    const G4Material*, const G4ParticleDefinition*);

  G4double TrueStepLength(G4double kinEnergy, G4double step,
			  const G4Material*, const G4ParticleDefinition* part);

  inline G4double EnergyAfterStep(G4double kinEnergy, G4double step, 
				  const G4Material*, 
                                  const G4String& particleName);

  inline G4double EnergyBeforeStep(G4double kinEnergy, G4double step, 
				   const G4Material*, 
                                   const G4String& particleName);

  G4double AverageScatteringAngle(G4double kinEnergy, G4double step, 
				  const G4Material*, 
				  const G4ParticleDefinition* part);

  inline G4double AverageScatteringAngle(G4double kinEnergy, G4double step, 
					 const G4Material*, 
					 const G4String& particleName);

  inline G4double ComputeTrueStep(const G4Material*, 
				  const G4ParticleDefinition* part, 
				  G4double kinEnergy, G4double stepLength);

  G4double EnergyDispersion(G4double kinEnergy, G4double step, 
			    const G4Material*, 
			    const G4ParticleDefinition*);

  inline G4double EnergyDispersion(G4double kinEnergy, G4double step, 
				   const G4Material*, 
                                   const G4String& particleName);

  inline void SetVerbose(G4int val);

  inline void SetMinKinEnergy(G4double);

  inline void SetMaxKinEnergy(G4double);

  inline void SetMaxEnergyTransfer(G4double);

  // hide assignment operator
  G4EnergyLossForExtrapolator & operator=
  (const G4EnergyLossForExtrapolator &right) = delete;
  G4EnergyLossForExtrapolator(const G4EnergyLossForExtrapolator&) = delete;
   
private:

  G4bool SetupKinematics(const G4ParticleDefinition*, const G4Material*, 
			 G4double kinEnergy);

  const G4ParticleDefinition* FindParticle(const G4String& name);

  inline G4double ComputeValue(G4double x, const G4PhysicsTable* table,
			       size_t idxMat);

  inline const G4PhysicsTable* GetPhysicsTable(ExtTableType type) const;

#ifdef G4MULTITHREADED
  static G4Mutex extrMutex;
#endif
  static G4TablesForExtrapolator* tables;

  const G4ParticleDefinition* currentParticle = nullptr;
  const G4ParticleDefinition* electron = nullptr;
  const G4ParticleDefinition* positron = nullptr;
  const G4ParticleDefinition* muonPlus = nullptr;
  const G4ParticleDefinition* muonMinus= nullptr;
  const G4ParticleDefinition* proton = nullptr;
  const G4Material* currentMaterial = nullptr;
   
  G4double electronDensity = 0.0;
  G4double radLength = 0.0;
  G4double charge2 = 0.0;
  G4double kineticEnergy = 0.0;
  G4double gam = 1.0;
  G4double bg2 = 0.0;
  G4double beta2 = 0.0;
  G4double tmax = 0.0;

  G4double linLossLimit = 0.01;
  G4double emin = 0.0;
  G4double emax = 0.0;
  G4double maxEnergyTransfer = 0.0;

  size_t index = 0;
  size_t  nmat = 0;
  G4int  nbins = 80;
  G4int  verbose = 0;

  G4bool isMaster = false;

  G4String currentParticleName = "";
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline const G4PhysicsTable* 
G4EnergyLossForExtrapolator::GetPhysicsTable(ExtTableType type) const
{
  return tables->GetPhysicsTable(type);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4EnergyLossForExtrapolator::EnergyAfterStep(G4double kinEnergy, 
					     G4double step, 
					     const G4Material* mat, 
					     const G4String& name)
{
  return EnergyAfterStep(kinEnergy,step,mat,FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4EnergyLossForExtrapolator::EnergyBeforeStep(G4double kinEnergy, 
					      G4double step, 
					      const G4Material* mat, 
					      const G4String& name)
{
  return EnergyBeforeStep(kinEnergy,step,mat,FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4EnergyLossForExtrapolator::AverageScatteringAngle(G4double kinEnergy, 
						    G4double step, 
						    const G4Material* mat, 
						    const G4String& name)
{
  return AverageScatteringAngle(kinEnergy,step,mat,FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4EnergyLossForExtrapolator::EnergyDispersion(G4double kinEnergy, 
					      G4double step, 
					      const G4Material* mat, 
					      const G4String& name)
{
  return EnergyDispersion(kinEnergy,step,mat,FindParticle(name));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4EnergyLossForExtrapolator::ComputeTrueStep(const G4Material* mat, 
					     const G4ParticleDefinition* part,
					     G4double kinEnergy, 
					     G4double stepLength)
{
  G4double theta = AverageScatteringAngle(kinEnergy,stepLength,mat,part);
  return stepLength*std::sqrt(1.0 + 0.625*theta*theta);
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4EnergyLossForExtrapolator::ComputeValue(G4double x, 
					  const G4PhysicsTable* table,
					  size_t idxMat)
{
  return (nullptr != table) ? ((*table)[idxMat])->Value(x, index) : 0.0; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolator::SetVerbose(G4int val) 
{
  verbose = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolator::SetMinKinEnergy(G4double val)
{
  emin = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolator::SetMaxKinEnergy(G4double val)
{
  emax = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4EnergyLossForExtrapolator::SetMaxEnergyTransfer(G4double val)
{
  maxEnergyTransfer = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

