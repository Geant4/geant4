//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4eeToHadrons.hh,v 1.1 2004/11/19 18:44:04 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eeToHadrons
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 12.08.2004
//
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
//
//
// Class Description:
//
// This class manages the process of e+ annihilation into hadrons
//

// -------------------------------------------------------------------
//

#ifndef G4eeToHadrons_h
#define G4eeToHadrons_h 1

#include "G4VEmProcess.hh"
#include "G4Positron.hh"
#include "G4eeToHadronsModel.hh"
#include <vector>

class G4eeCrossSections;

class G4eeToHadrons : public G4VEmProcess
{

public:

  G4eeToHadrons(const G4String& name = "ee2hadr");

  virtual ~G4eeToHadrons();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual std::vector<G4DynamicParticle*>* SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*);

  G4double CrossSection(G4double kineticEnergy, const G4MaterialCutsCouple* couple);
  // It returns the cross section of the process for energy/ material

  virtual G4double RecalculateLambda(G4double kinEnergy,
                               const G4MaterialCutsCouple* couple);

  virtual void PrintInfoDefinition();
  // Print out of the class parameters

  virtual G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);

  void SetCrossSecFactor(G4double fac);
  // Set the factor to artificially increase the crossSection (default 1)

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

  G4double GetMeanFreePath(const G4Track&,G4double,G4ForceCondition*);

  virtual void ResetNumberOfInteractionLengthLeft();

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dp);

private:

  G4double ComputeMeanFreePath(G4double kineticEnergy,
                         const G4MaterialCutsCouple* couple);

  std::vector<G4DynamicParticle*>* GenerateSecondaries(const G4DynamicParticle*);


  // hide assignment operator
  G4eeToHadrons & operator=(const G4eeToHadrons &right);
  G4eeToHadrons(const G4eeToHadrons&);

  G4eeCrossSections*               cross;

  std::vector<G4eeToHadronsModel*> models;
  std::vector<G4double>            ekinMin;
  std::vector<G4double>            ekinPeak;
  std::vector<G4double>            ekinMax;
  std::vector<G4double>            cumSum;

  G4double                         thKineticEnergy;
  G4double                         maxKineticEnergy;
  G4double                         mfpKineticEnergy;
  G4double                         preStepMFP;
  G4double                         preStepCS;
  G4double                         lambdaFactor;
  G4double                         csFactor;

  const G4MaterialCutsCouple*      currentCouple;
  
  G4int                            nModels;

  G4bool                           isInitialised;
  G4bool                           abovePeak;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4eeToHadrons::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadrons::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{
  return dp->GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4DynamicParticle*>* G4eeToHadrons::SecondariesPostStep(
                                                  G4VEmModel*,
                                            const G4MaterialCutsCouple*,
                                            const G4DynamicParticle* dp)
{
  std::vector<G4DynamicParticle*>* newp = 0;
  G4double kinEnergy = dp->GetKineticEnergy();
  if (kinEnergy > thKineticEnergy) newp = GenerateSecondaries(dp);
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadrons::GetMeanFreePath(const G4Track& track,
                                                     G4double,
                                                     G4ForceCondition* condition)
{
  *condition = NotForced;
  G4double kinEnergy = track.GetKineticEnergy();
  G4double x = DBL_MAX;
  if (kinEnergy > thKineticEnergy) 
    x = ComputeMeanFreePath(kinEnergy, track.GetMaterialCutsCouple());
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadrons::CrossSection(G4double kineticEnergy,
                                      const G4MaterialCutsCouple* couple)
{
  G4double cross = 0.0;
  if (kineticEnergy > thKineticEnergy) {
    for(G4int i=0; i<nModels; i++) {
      if(kineticEnergy >= ekinMin[i] && kineticEnergy <= ekinMax[i]) 
        cross += csFactor*(models[i])->CrossSection(couple,0,kineticEnergy,0.0,0.0);
      cumSum[i] = cross;
    } 
  }
  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eeToHadrons::RecalculateLambda(G4double kinEnergy, 
                                           const G4MaterialCutsCouple* couple)
{
  return CrossSection(kinEnergy, couple);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eeToHadrons::ResetNumberOfInteractionLengthLeft()
{
  currentCouple = 0;
  abovePeak     = false;
  preStepCS     = 0.0;
  G4VProcess::ResetNumberOfInteractionLengthLeft();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
