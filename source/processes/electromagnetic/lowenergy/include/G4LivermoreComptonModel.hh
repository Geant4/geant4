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
// $Id: G4LivermoreComptonModel.hh 82874 2014-07-15 15:25:29Z gcosmo $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Author: Alexander Bagulya
//         11 March 2012
//         on the base of G4LivermoreComptonModel
//
// History:
// --------

#ifndef G4LivermoreComptonModel_h
#define G4LivermoreComptonModel_h 1

#include "G4VEmModel.hh"
#include "G4LPhysicsFreeVector.hh"

class G4ParticleChangeForGamma;
class G4VAtomDeexcitation;
class G4ShellData;
class G4DopplerProfile;

class G4LivermoreComptonModel : public G4VEmModel
{
public:

  G4LivermoreComptonModel(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "LivermoreCompton");

  virtual ~G4LivermoreComptonModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual void InitialiseLocal(const G4ParticleDefinition*,
			       G4VEmModel* masterModel);

  virtual void InitialiseForElement(const G4ParticleDefinition*, G4int Z);

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

private:

  void ReadData(size_t Z, const char* path = 0);

  G4double ComputeScatteringFunction(G4double x, G4int Z);

  G4LivermoreComptonModel & operator=(const  G4LivermoreComptonModel &right);
  G4LivermoreComptonModel(const  G4LivermoreComptonModel&);

  G4bool isInitialised;
  G4int verboseLevel;
  
  G4ParticleChangeForGamma* fParticleChange;
  G4VAtomDeexcitation*      fAtomDeexcitation;

  static G4ShellData*       shellData;
  static G4DopplerProfile*  profileData;

  static G4int maxZ;
  static G4LPhysicsFreeVector* data[100];

  static const G4double ScatFuncFitParam[101][9];

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
