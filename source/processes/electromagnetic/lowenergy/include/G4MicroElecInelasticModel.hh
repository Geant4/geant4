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
//
// G4MicroElecInelasticModel.hh, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//
//          - Inelastic cross-sections of low energy electrons in silicon
//	    for the simulation of heavy ion tracks with theGeant4-DNA toolkit,
//	    NSS Conf. Record 2010, pp. 80-85
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	    NIM B, vol. 288, pp. 66-73, 2012.
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for protons and
//	    heavy ions in Si, NIM B, vol. 287, pp. 124-129, 2012.
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef G4MicroElecInelasticModel_h
#define G4MicroElecInelasticModel_h 1


#include "globals.hh"
#include "G4VEmModel.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4ProductionCutsTable.hh"

#include "G4MicroElecCrossSectionDataSet.hh"
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4GenericIon.hh"
#include "G4ParticleDefinition.hh"

#include "G4LogLogInterpolation.hh"

#include "G4MicroElecSiStructure.hh"
#include "G4VAtomDeexcitation.hh"
#include "G4NistManager.hh"

class G4MicroElecInelasticModel : public G4VEmModel
{

public:

  G4MicroElecInelasticModel(const G4ParticleDefinition* p = 0, 
		           const G4String& nam = "MicroElecInelasticModel");

  virtual ~G4MicroElecInelasticModel();

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

  double DifferentialCrossSection(G4ParticleDefinition * aParticleDefinition, G4double k, G4double energyTransfer, G4int shell);

  G4double TransferedEnergy(G4ParticleDefinition * aParticleDefinition,
                            G4double incomingParticleEnergy, G4int shell, G4double random) ;
  
  inline void SelectFasterComputation(G4bool input); 
  
protected:

  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
  
  G4bool fasterCode;

  //deexcitation manager to produce fluo photns and e-
  G4VAtomDeexcitation*      fAtomDeexcitation;

  G4Material* nistSi;

  std::map<G4String,G4double,std::less<G4String> > lowEnergyLimit;
  std::map<G4String,G4double,std::less<G4String> > highEnergyLimit;

  G4bool isInitialised;
  G4int verboseLevel;
  
  // Cross section

  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;

  typedef std::map<G4String,G4MicroElecCrossSectionDataSet*,std::less<G4String> > MapData;
  MapData tableData;
  
  // Final state
  
  G4MicroElecSiStructure SiStructure;

  G4double RandomizeEjectedElectronEnergy(G4ParticleDefinition * aParticleDefinition, G4double incomingParticleEnergy, G4int shell) ;

  G4double RandomizeEjectedElectronEnergyFromCumulatedDcs(G4ParticleDefinition * aParticleDefinition, G4double incomingParticleEnergy, G4int shell) ;

  G4double Interpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
   
  G4double QuadInterpolator( G4double e11, 
			     G4double e12, 
			     G4double e21, 
			     G4double e22, 
			     G4double x11,
			     G4double x12, 
			     G4double x21, 
			     G4double x22, 
			     G4double t1, 
			     G4double t2, 
			     G4double t, 
			     G4double e);

  typedef std::map<double, std::map<double, double> > TriDimensionMap;
  
  TriDimensionMap eDiffCrossSectionData[7];
  TriDimensionMap eNrjTransfData[7]; // for cumulated dcs
  
  TriDimensionMap pDiffCrossSectionData[7];
  TriDimensionMap pNrjTransfData[7]; // for cumulated dcs
  
  std::vector<double> eTdummyVec;
  std::vector<double> pTdummyVec;

  typedef std::map<double, std::vector<double> > VecMap;
  
  VecMap eVecm;
  VecMap pVecm;
  
  VecMap eProbaShellMap[7]; // for cumulated dcs
  VecMap pProbaShellMap[7]; // for cumulated dcs
  
  // Partial cross section
  
  G4int RandomSelect(G4double energy,const G4String& particle );
   
  //
   
  G4MicroElecInelasticModel & operator=(const  G4MicroElecInelasticModel &right);
  G4MicroElecInelasticModel(const  G4MicroElecInelasticModel&);

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4MicroElecInelasticModel::SelectFasterComputation (G4bool input)
{ 
    fasterCode = input; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
