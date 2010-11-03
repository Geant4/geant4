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
// $Id: G4DNARuddIonisationExtendedModel.hh,v 1.1 2010-11-03 10:44:26 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4DNARuddIonisationExtendedModel_h
#define G4DNARuddIonisationExtendedModel_h 1

#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4DNAGenericIonsManager.hh"
#include "G4DNACrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4LogLogInterpolation.hh"

#include "G4WaterIonisationStructure.hh"

class G4DNARuddIonisationExtendedModel : public G4VEmModel
{

public:

  G4DNARuddIonisationExtendedModel(const G4ParticleDefinition* p = 0, 
		           const G4String& nam = "DNARuddIonisationExtendedModel");

  virtual ~G4DNARuddIonisationExtendedModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double CrossSectionPerVolume(  const G4Material* material,
					   const G4ParticleDefinition* p,
					   G4double ekin,
					   G4double emin,
					   G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:

  std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
  std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;
// ZF 26-10-2010
  std::map<G4double,G4double> lowEnergyLimitForA, lowEnergyLimitOfModelForA, killBelowEnergyForA;

  G4bool isInitialised;
  G4int verboseLevel;
  
  // Cross section

  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;

  typedef std::map<G4String,G4DNACrossSectionDataSet*,std::less<G4String> > MapData;
  MapData tableData;
  
  // Final state
  
  G4WaterIonisationStructure waterStructure;

  G4double RandomizeEjectedElectronEnergy(G4ParticleDefinition* particleDefinition, 
					  G4double incomingParticleEnergy, 
					  G4int shell);

  void RandomizeEjectedElectronDirection(G4ParticleDefinition* particleDefinition, 
					 G4double incomingParticleEnergy, 
					 G4double outgoingParticleEnergy, 
					 G4double & cosTheta, 
					 G4double & phi, G4int shell);

  G4double RejectionFunction(G4ParticleDefinition* particle, 
							      G4double k, 
							      G4double proposed_ws, 
							      G4int ionizationLevelIndex);

  G4double ProposedSampledEnergy(G4ParticleDefinition* particle, 
							      G4double k,  
							      G4int ionizationLevelIndex);
   
  G4double CorrectionFactor(G4ParticleDefinition* particleDefinition, G4double k, G4int shell);

  G4double S_1s(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveChg, 
		G4double shellNumber);

  G4double S_2s(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveChg, 
		G4double shellNumber);


  G4double S_2p(G4double t, 
		G4double energyTransferred, 
		G4double slaterEffectiveChg, 
		G4double shellNumber);

  G4double R(G4double t, 
	     G4double energyTransferred, 
	     G4double slaterEffectiveChg, 
	     G4double shellNumber) ;

  G4double slaterEffectiveCharge[3];
  G4double sCoefficient[3];
  
  // Partial cross section
  
  G4double PartialCrossSection(const G4Track& track);
			
  G4double Sum(G4double energy, const G4String& particle);

  G4int RandomSelect(G4double energy,const G4String& particle );
  
  // Test water material 
   
  G4bool flagMaterialIsWater;
  G4double densityWater;
   
  //
   
  G4DNARuddIonisationExtendedModel & operator=(const  G4DNARuddIonisationExtendedModel &right);
  G4DNARuddIonisationExtendedModel(const  G4DNARuddIonisationExtendedModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
