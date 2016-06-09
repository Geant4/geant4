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
// $Id: G4eplusAnnihilation70.hh,v 1.3 2004/11/10 08:53:19 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-01 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eplusAnnihilation70
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
//
//
// Class Description:
//
// This class manages the process of e+ annihilation into 2 gammas
//

// -------------------------------------------------------------------
//

#ifndef G4eplusAnnihilation70_h
#define G4eplusAnnihilation70_h 1

#include "G4VEmProcess.hh"
#include "G4Positron.hh"

class G4eplusAnnihilation70 : public G4VEmProcess
{

public:

  G4eplusAnnihilation70(const G4String& name = "annihil");

  virtual ~G4eplusAnnihilation70();

  G4bool IsApplicable(const G4ParticleDefinition& p);

  virtual std::vector<G4DynamicParticle*>* SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*);

  virtual G4VParticleChange* AtRestDoIt(
                             const G4Track& track,
                             const G4Step& stepData);

  virtual G4double AtRestGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4ForceCondition* condition
                            );

  void PrintInfoDefinition();
  // Print out of the class parameters

  G4PhysicsVector* LambdaPhysicsVector(const G4MaterialCutsCouple*);


protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

  virtual G4double MaxSecondaryEnergy(const G4DynamicParticle* dp);

private:

 // hide assignment operator
  G4eplusAnnihilation70 & operator=(const G4eplusAnnihilation70 &right);
  G4eplusAnnihilation70(const G4eplusAnnihilation70&);
  
  G4bool isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4eplusAnnihilation70::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eplusAnnihilation70::MaxSecondaryEnergy(const G4DynamicParticle* dp)
{
  return dp->GetKineticEnergy();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eplusAnnihilation70::AtRestGetPhysicalInteractionLength(
                              const G4Track&, G4ForceCondition* condition)
{
  *condition = NotForced;
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmModel.hh"

inline std::vector<G4DynamicParticle*>* G4eplusAnnihilation70::SecondariesPostStep(
                                                  G4VEmModel* model,
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* dp)
{
  fParticleChange.ProposeTrackStatus(fStopAndKill);
  return model->SampleSecondaries(couple, dp, 0.0, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
