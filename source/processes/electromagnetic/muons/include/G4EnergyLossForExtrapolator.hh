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
// $Id: G4EnergyLossForExtrapolator.hh 97392 2016-06-02 10:10:32Z gcosmo $
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

class G4ParticleDefinition;
class G4Material;
class G4MaterialCutsCouple;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4EnergyLossForExtrapolator 
{
public:

  explicit G4EnergyLossForExtrapolator(G4int verb = 1);

  ~G4EnergyLossForExtrapolator();

  G4double ComputeDEDX(G4double kinEnergy, const G4ParticleDefinition*);

  G4double ComputeRange(G4double kinEnergy, const G4ParticleDefinition*);

  G4double ComputeEnergy(G4double range, const G4ParticleDefinition*);

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

  inline G4double AverageScatteringAngle(G4double kinEnergy, G4double step, 
					 const G4Material*, 
					 const G4ParticleDefinition* part);

  inline G4double AverageScatteringAngle(G4double kinEnergy, G4double step, 
					 const G4Material*, 
					 const G4String& particleName);

  inline G4double ComputeTrueStep(const G4Material*, 
				  const G4ParticleDefinition* part, 
				  G4double kinEnergy, G4double stepLength);

  inline G4double EnergyDispersion(G4double kinEnergy, G4double step, 
				   const G4Material*, 
                                   const G4ParticleDefinition*);

  inline G4double EnergyDispersion(G4double kinEnergy, G4double step, 
				   const G4Material*, 
                                   const G4String& particleName);

  inline void SetVerbose(G4int val);

  inline void SetMinKinEnergy(G4double);

  inline void SetMaxKinEnergy(G4double);

  inline void SetMaxEnergyTransfer(G4double);
   
private:

  void Initialisation();

  void BuildTables();

  G4bool SetupKinematics(const G4ParticleDefinition*, const G4Material*, 
			 G4double kinEnergy);

  const G4ParticleDefinition* FindParticle(const G4String& name);

  inline G4double ComputeValue(G4double x, const G4PhysicsTable* table,
			       size_t idx);

  inline const G4PhysicsTable* GetPhysicsTable(ExtTableType type) const;

  // hide assignment operator
  G4EnergyLossForExtrapolator & operator=(const G4EnergyLossForExtrapolator &right);
  G4EnergyLossForExtrapolator(const G4EnergyLossForExtrapolator&);

  static G4TablesForExtrapolator* tables;

  const G4ParticleDefinition* currentParticle;
  const G4ParticleDefinition* electron;
  const G4ParticleDefinition* positron;
  const G4ParticleDefinition* muonPlus;
  const G4ParticleDefinition* muonMinus;
  const G4ParticleDefinition* proton;

  G4String currentParticleName;
  
  size_t idxDedxElectron;
  size_t idxDedxPositron;
  size_t idxDedxMuon;
  size_t idxDedxProton;
  size_t idxRangeElectron;
  size_t idxRangePositron;
  size_t idxRangeMuon;
  size_t idxRangeProton;
  size_t idxInvRangeElectron;
  size_t idxInvRangePositron;
  size_t idxInvRangeMuon;
  size_t idxInvRangeProton;
  size_t idxMscElectron;
 
  const G4Material* currentMaterial;
  G4int       index;

  G4double    electronDensity;
  G4double    radLength;
  G4double    mass;
  G4double    charge2;
  G4double    kineticEnergy;
  G4double    gam;
  G4double    bg2;
  G4double    beta2;
  G4double    tmax;

  G4double    linLossLimit;
  G4double    emin;
  G4double    emax;
  G4double    maxEnergyTransfer;

  G4int       nbins;
  G4int       nmat;
  G4int       verbose;
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

inline G4double G4EnergyLossForExtrapolator::AverageScatteringAngle(
                        G4double kinEnergy, 
			G4double stepLength, 
			const G4Material* mat, 
			const G4ParticleDefinition* part)
{
  G4double theta = 0.0;
  if(SetupKinematics(part, mat, kinEnergy)) {
    G4double t = stepLength/radLength;
    G4double y = std::max(0.001, t); 
    theta = 19.23*CLHEP::MeV*std::sqrt(charge2*t)*(1.0 + 0.038*G4Log(y))
      /(beta2*gam*mass);
  }
  return theta;
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
G4EnergyLossForExtrapolator::EnergyDispersion(G4double kinEnergy, 
					      G4double stepLength, 
					      const G4Material* mat, 
					      const G4ParticleDefinition* part)
{
  G4double sig2 = 0.0;
  if(SetupKinematics(part, mat, kinEnergy)) {
    G4double step = ComputeTrueStep(mat,part,kinEnergy,stepLength);
    sig2 = (1.0/beta2 - 0.5)
      *CLHEP::twopi_mc2_rcl2*tmax*step*electronDensity*charge2;
  }
  return sig2;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double 
G4EnergyLossForExtrapolator::ComputeValue(G4double x, 
					  const G4PhysicsTable* table,
					  size_t idx)
{
  G4double res = 0.0;
  if(table) { res = ((*table)[index])->Value(x, idx); }
  return res;
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

