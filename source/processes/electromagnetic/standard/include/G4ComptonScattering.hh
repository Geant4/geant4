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
//
// $Id: G4ComptonScattering.hh,v 1.15 2005/05/12 11:06:42 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//------------------ G4ComptonScattering physics process -----------------------
//                   by Michel Maire, April 1996
//
// 10-06-96, updated by M.Maire
// 21-06-96, SetCuts implementation, M.Maire
// 06-01-97, crossection table + meanfreepath table, M.Maire
// 17-02-97, New Physics scheme
// 25-02-97, GetMeanFreePath() now is public function
// 12-03-97, new physics scheme again
// 13-08-98, new methods SetBining()  PrintInfo()
// 03-08-01, new methods Store/Retrieve PhysicsTable (mma)
// 06-08-01, BuildThePhysicsTable() called from constructor (mma)
// 19-09-01, come back to previous ProcessName "compt"
// 20-09-01, DoIt: fminimalEnergy = 1*eV (mma)
// 01-10-01, come back to BuildPhysicsTable(const G4ParticleDefinition&)
// 13-08-04, suppress icc file; make public ComputeCrossSectionPerAtom()   (mma)
// 09-11-04, Remove Retrieve tables (V.Ivantchenko)
// 15-03-05, Redesign - use G4VEmProcess interface (V.Ivantchenko)
// 04-05-05, Make class to be default (V.Ivanchenko)
// -----------------------------------------------------------------------------

// class description
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4ComptonScattering_h
#define G4ComptonScattering_h 1

#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4Gamma.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4ParticleDefinition;
class G4VEmModel;
class G4MaterialCutsCouple;
class G4DynamicParticle;

class G4ComptonScattering : public G4VEmProcess

{
public:  // with description

  G4ComptonScattering(const G4String& processName ="compt",
		      G4ProcessType type = fElectromagnetic);

  virtual ~G4ComptonScattering();

  // true for Gamma only.  
  G4bool IsApplicable(const G4ParticleDefinition&);

  // Print few lines of informations about the process: validity range,
  virtual void PrintInfo();

  void SetModel(const G4String& name);

protected:

  virtual void InitialiseProcess(const G4ParticleDefinition*);

  std::vector<G4DynamicParticle*>* SecondariesPostStep(
                                   G4VEmModel*,
                             const G4MaterialCutsCouple*,
                             const G4DynamicParticle*);

private:
  
  // hide assignment operator as private 
  G4ComptonScattering& operator=(const G4ComptonScattering &right);
  G4ComptonScattering(const G4ComptonScattering& );
     
  G4bool          isInitialised;
  G4VEmModel*     selectedModel;
  G4int           mType;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4bool G4ComptonScattering::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4Gamma::Gamma());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline std::vector<G4DynamicParticle*>* G4ComptonScattering::SecondariesPostStep(
                                   G4VEmModel* model,
                             const G4MaterialCutsCouple* couple,
                             const G4DynamicParticle* dp)
{ 
  return model->SampleSecondaries(couple, dp);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
#endif
 
