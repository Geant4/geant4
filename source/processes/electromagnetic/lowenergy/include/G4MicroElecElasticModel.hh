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
// G4MicroElecElasticModel.hh, 2011/08/29 A.Valentin, M. Raine
//
// Based on the following publications
//	    - Geant4 physics processes for microdosimetry simulation:
//	    very low energy electromagnetic models for electrons in Si,
//	     NIM B, vol. 288, pp. 66 - 73, 2012.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef G4MicroElecElasticModel_h
#define G4MicroElecElasticModel_h 1

#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4MicroElecCrossSectionDataSet.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"

class G4MicroElecElasticModel : public G4VEmModel
{

public:
  G4MicroElecElasticModel(const G4ParticleDefinition* p = nullptr, 
		          const G4String& nam = "MicroElecElasticModel");
  virtual ~G4MicroElecElasticModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double CrossSectionPerVolume(const G4Material* material,
				 const G4ParticleDefinition* p,
				 G4double ekin,
				 G4double emin,
				 G4double emax) override;
  
  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  inline void SetKillBelowThreshold (G4double threshold);		 
  G4double GetKillBelowThreshold () { return killBelowEnergy; }	

  G4MicroElecElasticModel & operator=(const  G4MicroElecElasticModel &right) = delete;
  G4MicroElecElasticModel(const  G4MicroElecElasticModel&) = delete;

protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma;

private:
  // Final state
  G4double Theta(G4ParticleDefinition * aParticleDefinition, G4double k, G4double integrDiff);
  G4double LinLinInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
  G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);   
  G4double LinLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);   
  G4double QuadInterpolator(G4double e11, 
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
  G4double RandomizeCosTheta(G4double k);

  // Cross section 
  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;

  typedef std::map<G4String,G4MicroElecCrossSectionDataSet*,std::less<G4String> > MapData;
  MapData tableData;
  
  typedef std::map<G4double, std::map<G4double, G4double> > TriDimensionMap;
  TriDimensionMap eDiffCrossSectionData;
  std::vector<G4double> eTdummyVec;

  typedef std::map<G4double, std::vector<G4double> > VecMap;
  VecMap eVecm;

  G4Material* nistSi;
  G4double killBelowEnergy;  
  G4double lowEnergyLimit;  
  G4double lowEnergyLimitOfModel;  
  G4double highEnergyLimit; 
  G4int verboseLevel;
  G4bool isInitialised;    
};

inline void G4MicroElecElasticModel::SetKillBelowThreshold (G4double threshold) 
{ 
    killBelowEnergy = threshold; 
    
    if (threshold < 5*CLHEP::eV)
    {
       G4Exception ("*** WARNING : the G4MicroElecElasticModel class is not validated below 5 eV !","",JustWarning,"") ;
       threshold = 5*CLHEP::eV;
    }
             
}		 


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
