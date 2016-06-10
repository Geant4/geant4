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
// $Id: G4LivermoreComptonModifiedModel.hh 82874 2014-07-15 15:25:29Z gcosmo $
//
// Author: Sebastien Incerti
//         30 October 2008
//         on base of G4LowEnergyCompton developed by A.Forti and M.G.Pia
//
// Modified:
// --------
// 30 May 2011   V Ivanchenko Migration to model design for deexcitation

#ifndef G4LivermoreComptonModifiedModel_h
#define G4LivermoreComptonModifiedModel_h 1

#include "G4VEmModel.hh"
#include "G4ShellData.hh"
#include "G4DopplerProfile.hh"

class G4ParticleChangeForGamma;
class G4VCrossSectionHandler;
class G4VAtomDeexcitation;
class G4VEMDataSet;

class G4LivermoreComptonModifiedModel : public G4VEmModel
{

public:

  G4LivermoreComptonModifiedModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "LivermoreModifiedCompton");

  virtual ~G4LivermoreComptonModifiedModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom( const G4ParticleDefinition*,
                                               G4double kinEnergy, 
                                               G4double Z, 
                                               G4double A=0, 
                                               G4double cut=0,
                                               G4double emax=DBL_MAX );

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

protected:

  G4ParticleChangeForGamma* fParticleChange;

private:
  G4bool isInitialised;
  G4int verboseLevel;
  
  G4VEMDataSet* scatterFunctionData;
  G4VCrossSectionHandler* crossSectionHandler;

  G4VAtomDeexcitation*    fAtomDeexcitation;

  G4ShellData shellData;
  G4DopplerProfile profileData;

  G4LivermoreComptonModifiedModel & operator=(const  G4LivermoreComptonModifiedModel &right);
  G4LivermoreComptonModifiedModel(const  G4LivermoreComptonModifiedModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
