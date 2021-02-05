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
// G4MicroElecLOPhononModel.hh, 
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
//
//	- Q.Gibaru, C.Inguimbert, P.Caron, M.Raine, D.Lambert, J.Puech, 
//	      Geant4 physics processes for microdosimetry and secondary electron emission simulation : 
//	      Extension of MicroElec to very low energies and new materials
//	      NIM B, 2020, in review.
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef G4MICROELECLOPHONONMODEL_HH 
#define G4MICROELECLOPHONONMODEL_HH 1 

#include "G4Step.hh"
#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Material.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleChangeForGamma.hh"

class G4MicroElecLOPhononModel : public G4VEmModel
{
public:
  G4MicroElecLOPhononModel(const G4ParticleDefinition*p = nullptr,
		  const G4String& nam = "G4MicroElecLOPhononModel");
  ~G4MicroElecLOPhononModel() override;
  
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
  
protected:
  G4ParticleChangeForGamma* fParticleChangeForGamma;
  
private:
  G4MicroElecLOPhononModel & operator=(const  G4MicroElecLOPhononModel &right);
  G4MicroElecLOPhononModel(const  G4MicroElecLOPhononModel&);
  
  G4bool Interband = false;
  G4bool isInitialised = false;
  G4bool absor = false;
  G4double phononEnergy = 0.;
};

#endif 
