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
// G4MicroElecElasticModel_new.hh, 2011/08/29 A.Valentin, M. Raine are with CEA [a]
//                   	    	2020/05/20 P. Caron, C. Inguimbert are with ONERA [b] 
//				       	   Q. Gibaru is with CEA [a], ONERA [b] and CNES [c]
//				           M. Raine and D. Lambert are with CEA [a]
//
// A part of this work has been funded by the French space agency(CNES[c])
// [a] CEA, DAM, DIF - 91297 ARPAJON, France
// [b] ONERA - DPHY, 2 avenue E.Belin, 31055 Toulouse, France
// [c] CNES, 18 av.E.Belin, 31401 Toulouse CEDEX, France
//
// Based on the following publications
//	- A.Valentin, M. Raine, 
//		Inelastic cross-sections of low energy electrons in silicon
//	      for the simulation of heavy ion tracks with the Geant4-DNA toolkit,
//	      NSS Conf. Record 2010, pp. 80-85
//             https://doi.org/10.1109/NSSMIC.2010.5873720
//
//      - A.Valentin, M. Raine, M.Gaillardin, P.Paillet
//	      Geant4 physics processes for microdosimetry simulation:
//	      very low energy electromagnetic models for electrons in Silicon,
//             https://doi.org/10.1016/j.nimb.2012.06.007
//	      NIM B, vol. 288, pp. 66-73, 2012, part A
//	      heavy ions in Si, NIM B, vol. 287, pp. 124-129, 2012, part B
//             https://doi.org/10.1016/j.nimb.2012.07.028
//
//	- M. Raine, M. Gaillardin, P. Paillet
//	      Geant4 physics processes for silicon microdosimetry simulation: 
//	      Improvements and extension of the energy-range validity up to 10 GeV/nucleon
//	      NIM B, vol. 325, pp. 97-100, 2014
//             https://doi.org/10.1016/j.nimb.2014.01.014
//
//      - J. Pierron, C. Inguimbert, M. Belhaj, T. Gineste, J. Puech, M. Raine
//	      Electron emission yield for low energy electrons: 
//	      Monte Carlo simulation and experimental comparison for Al, Ag, and Si
//	      Journal of Applied Physics 121 (2017) 215107. 
//               https://doi.org/10.1063/1.4984761
//
//      - P. Caron,
//	      Study of Electron-Induced Single-Event Upset in Integrated Memory Devices
//	      PHD, 16th October 2019
//
//	- Q.Gibaru, C.Inguimbert, P.Caron, M.Raine, D.Lambert, J.Puech, 
//	      Geant4 physics processes for microdosimetry and secondary electron emission simulation : 
//	      Extension of MicroElec to very low energies and new materials
//	      NIM B, 2020, in review.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef G4MICROELECELASTICMODEL_NEW_HH
#define G4MICROELECELASTICMODEL_NEW_HH 1

#include <map>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4MicroElecMaterialStructure.hh"
#include "G4MicroElecCrossSectionDataSet_new.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4ParticleChangeForGamma.hh"
#include "G4LogLogInterpolation.hh"
#include "G4ProductionCutsTable.hh"
#include "G4NistManager.hh"

class G4MicroElecElasticModel_new : public G4VEmModel
{

public:
  G4MicroElecElasticModel_new(const G4ParticleDefinition* p = 0, 
		          const G4String& nam = "MicroElecElasticModel");
  ~G4MicroElecElasticModel_new() override;
  
  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;

  G4double CrossSectionPerVolume(const G4Material* material,
				 const G4ParticleDefinition* p,
				 G4double ekin,
				 G4double emin,
				 G4double emax) override;

  G4double AcousticCrossSectionPerVolume(G4double ekin, G4double kbz, G4double rho,
					 G4double cs, G4double Aac, G4double Eac,
					 G4double prefactor);

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double tmin,
			 G4double maxEnergy) override;

  void SetKillBelowThreshold (G4double threshold);		 

  G4double GetKillBelowThreshold () { return killBelowEnergy; }	

  G4double DamageEnergy(G4double T,G4double A, G4double Z);

protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma;
  
private:

  G4MicroElecElasticModel_new & operator=(const  G4MicroElecElasticModel_new &right);
  G4MicroElecElasticModel_new(const  G4MicroElecElasticModel_new&);

  // Final state
  G4double Theta(G4ParticleDefinition * aParticleDefinition, G4double k, G4double integrDiff);
  G4double LinLinInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
  G4double LogLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
  G4double LinLogInterpolate(G4double e1, G4double e2, G4double e, G4double xs1, G4double xs2);
  G4double QuadInterpolator(G4double e11, G4double e12, G4double e21, G4double e22, 
			    G4double x11, G4double x12, G4double x21, G4double x22, 
			    G4double t1, G4double t2, G4double t, G4double e);

  G4double RandomizeCosTheta(G4double k);

  G4Material* nistSi = nullptr;
  G4double killBelowEnergy;  
  G4double lowEnergyLimit;  
  G4double lowEnergyLimitOfModel;  
  G4double highEnergyLimit; 
  G4bool isInitialised;
  G4int verboseLevel;
  // Cross section
  typedef std::map<G4String,G4String,std::less<G4String> > MapFile;
  MapFile tableFile;
  typedef std::map<G4String,G4MicroElecCrossSectionDataSet_new*,std::less<G4String> > MapData;
  //MapData tableData;
  
  typedef std::map<G4String, MapData*, std::less<G4String> > TCSMap;
  TCSMap tableTCS;

  //Maps for multilayers
  typedef std::map<G4double, std::map<G4double, G4double> > TriDimensionMap;

  typedef std::map<G4String, TriDimensionMap* > ThetaMap;
  ThetaMap thetaDataStorage; //Storage of angles (cumulated)

  typedef std::map<G4String, std::vector<G4double>* > energyMap;
  energyMap eIncidentEnergyStorage;

  typedef std::map<G4double, std::vector<G4double> > VecMap;

  typedef std::map<G4String, VecMap* > ProbaMap;
  ProbaMap eProbaStorage; //Storage of probabilities for cumulated sections

  typedef std::map<G4String, G4MicroElecMaterialStructure*, std::less<G4String> > MapStructure;

  MapStructure tableMaterialsStructures; //Structures of all materials simulated

  G4MicroElecMaterialStructure* currentMaterialStructure = nullptr;
  typedef std::map<G4String, G4double, std::less<G4String> > MapEnergy;
  MapEnergy lowEnergyLimitTable;
  MapEnergy highEnergyLimitTable;
  MapEnergy workFunctionTable;
      
  G4bool killElectron, acousticModelEnabled;
  G4String currentMaterialName;
  G4bool isOkToBeInitialised;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
