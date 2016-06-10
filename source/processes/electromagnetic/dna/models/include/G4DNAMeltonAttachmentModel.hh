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
// $Id: G4DNAMeltonAttachmentModel.hh 70171 2013-05-24 13:34:18Z gcosmo $
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Created by Z. Francis

#ifndef G4DNAMeltonAttachmentModel_h
#define G4DNAMeltonAttachmentModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4DNAGenericIonsManager.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4LogLogInterpolation.hh"
#include "G4NistManager.hh"

class G4DNAMeltonAttachmentModel : public G4VEmModel
{

public:

  G4DNAMeltonAttachmentModel(const G4ParticleDefinition* p = 0, 
		             const G4String& nam = "DNAMeltonAttachmentModel");

  virtual ~G4DNAMeltonAttachmentModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(const G4Material* material,
					 const G4ParticleDefinition* p,
					 G4double ekin,
					 G4double emin,
					 G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  inline void SetDissociationFlag(G4bool);
  inline G4bool GetDissociationFlag();

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
  // Water density table
  const std::vector<G4double>* fpWaterDensity;

  G4double lowEnergyLimit;
  G4double highEnergyLimit;
  G4double lowEnergyLimitOfModel;
  G4bool isInitialised;
  G4int verboseLevel;
  G4bool fDissociationFlag;
  
  // Cross section
 
  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;

  typedef std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> > MapData;
  MapData tableData;
  
  //
  
  G4DNAMeltonAttachmentModel & operator=(const  G4DNAMeltonAttachmentModel &right);
  G4DNAMeltonAttachmentModel(const  G4DNAMeltonAttachmentModel&);

};

inline void G4DNAMeltonAttachmentModel::SetDissociationFlag(G4bool flag)
{
    fDissociationFlag = flag;
}

inline G4bool G4DNAMeltonAttachmentModel::GetDissociationFlag()
{
    return fDissociationFlag;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
