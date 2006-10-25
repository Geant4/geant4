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
// $Id: G4mplIonisation.hh,v 1.1 2006-10-25 17:37:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4mplIonisation
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 25.08.2005
//
// Modifications:
//
//
// Class Description:
//
// This class manages the ionisation process for a magnetic monopole
// it inherites from G4VContinuousDiscreteProcess via G4VEnergyLossProcess.
//

// -------------------------------------------------------------------
//

#ifndef G4mplIonisation_h
#define G4mplIonisation_h 1

#include "G4VEnergyLossProcess.hh"
#include "globals.hh"
#include "G4VEmModel.hh"

class G4Material;
class G4VEmFluctuationModel;

class G4mplIonisation : public G4VEnergyLossProcess
{

public:

  G4mplIonisation(const G4String& name = "mplIoni");

  virtual ~G4mplIonisation();

  G4bool IsApplicable(const G4ParticleDefinition& p);

  // Print out of the class parameters
  virtual void PrintInfo();

protected:

  std::vector<G4DynamicParticle*>*  SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*,
                                   G4double&);

  virtual void InitialiseEnergyLossProcess(const G4ParticleDefinition*,
					   const G4ParticleDefinition*);

private:

  // hide assignment operator
  G4mplIonisation & operator=(const G4mplIonisation &right);
  G4mplIonisation(const G4mplIonisation&);

  G4bool                      isInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4mplIonisation::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetParticleName() == "monopole");
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4DynamicParticle*>* G4mplIonisation::SecondariesPostStep(
                                                  G4VEmModel*,
                                            const G4MaterialCutsCouple*,
                                            const G4DynamicParticle*,
                                                  G4double&)
{
  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
