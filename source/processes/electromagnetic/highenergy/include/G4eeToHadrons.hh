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
// $Id: G4eeToHadrons.hh,v 1.4 2005/05/18 10:12:32 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
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
#include "G4eeToHadronsMultiModel.hh"

class G4eeToHadrons : public G4VEmProcess
{

public:

  G4eeToHadrons(const G4String& name = "ee2hadr");

  virtual ~G4eeToHadrons();

  virtual G4bool IsApplicable(const G4ParticleDefinition& p);

  // Print out of the class parameters
  void PrintInfo();

  // Set the factor to artificially increase the crossSection (default 1)
  void SetCrossSecFactor(G4double fac);

protected:

  void InitialiseProcess(const G4ParticleDefinition*);

  std::vector<G4DynamicParticle*>* SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*);

private:

  std::vector<G4DynamicParticle*>* GenerateSecondaries(const G4DynamicParticle*);

  // hide assignment operator
  G4eeToHadrons & operator=(const G4eeToHadrons &right);
  G4eeToHadrons(const G4eeToHadrons&);

  G4eeToHadronsMultiModel*  multimodel;
  G4bool                    isInitialised;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4eeToHadrons::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Positron::Positron());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4DynamicParticle*>* G4eeToHadrons::SecondariesPostStep(
                                                  G4VEmModel*,
                                            const G4MaterialCutsCouple* couple,
                                            const G4DynamicParticle* dp)
{
  std::vector<G4DynamicParticle*>* newp = multimodel->SampleSecondaries(couple, dp);
  if(newp) fParticleChange.ProposeTrackStatus(fStopAndKill);
  return newp;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
