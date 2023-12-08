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
// G4MicroElecSurface.hh, 
//                   	    	2020/05/20 P. Caron, C. Inguimbert are with ONERA [b] 
//				       	   Q. Gibaru is with CEA [a], ONERA [b] and CNES [c]
//				           D. Lambert is with CEA [a]
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
// Based on:
//		-the class G4OpBoundaryProcess.cc for the surface crossing of 
//		optical photons. 
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
 
#ifndef G4MicroElecSurface_h 
#define G4MicroElecSurface_h 1 
 
///////////// 
// Includes 
///////////// 
 
#include "globals.hh" 
#include "templates.hh" 
#include "geomdefs.hh" 
#include "Randomize.hh" 
#include "G4ProductionCutsTable.hh" 
#include "G4RandomTools.hh" 
#include "G4RandomDirection.hh" 
#include "G4MicroElecMaterialStructure.hh" 
#include "G4Step.hh" 
#include "G4VDiscreteProcess.hh" 
#include "G4DynamicParticle.hh" 
#include "G4Material.hh" 
#include "G4LogicalBorderSurface.hh" 
#include "G4LogicalSkinSurface.hh" 
#include "G4OpticalPhoton.hh"
#include "G4Electron.hh"
#include "G4TransportationManager.hh" 
  
// Class Description: 
// Discrete Process -- reflection/refraction at interfaces for electrons. 
// Class inherits publicly from G4VDiscreteProcess. 
// Class Description - End: 
 
///////////////////// 
// Class Definition 
///////////////////// 
 
enum G4MicroElecSurfaceStatus {  UndefinedSurf, 
				 NotAtBoundarySurf, 
				 SameMaterialSurf,
				 StepTooSmallSurf }; 
 
class G4MicroElecSurface : public G4VDiscreteProcess 
{ 
 public: 
  explicit G4MicroElecSurface(const G4String& processName = "MicroElecSurface", 
			      G4ProcessType type = fElectromagnetic);
  
  ~G4MicroElecSurface() override; 
 
  G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override; 
  // Returns true -> 'is applicable' only for an electron. 
  
  void SetFlagFranchissement(); 

  G4double GetMeanFreePath(const G4Track& , 
			   G4double , 
			   G4ForceCondition* condition) override; 
  // Returns infinity; i. e. the process does not limit the step, 
  // but sets the 'Forced' condition for the DoIt to be invoked at 
  // every step. However, only at a boundary will any action be 
  // taken. 
    
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, 
				  const G4Step&  aStep) override; 
  // This is the method implementing boundary processes. 

  void BuildPhysicsTable(const G4ParticleDefinition&) override;
  // Initialisation
 
  G4MicroElecSurfaceStatus GetStatus() const; 
  // Returns the current status. 

  G4MicroElecSurface(const G4MicroElecSurface &right) = delete; 
  G4MicroElecSurface& operator=(const G4MicroElecSurface &right) = delete; 
    
  void Initialise();

private:
   // Returns the incident angle of electron	
  G4double GetIncidentAngle(); 
  G4ThreeVector Reflexion(const G4StepPoint* PostStepPoint); 
  
  // private elements
  typedef std::map<G4String, G4double, std::less<G4String> > WorkFunctionTable;
  WorkFunctionTable tableWF; //Table of all materials simulated 
 
  G4double theParticleMomentum;
  G4ThreeVector oldMomentum, previousMomentum; 
  G4ThreeVector theGlobalNormal; 
  G4ThreeVector theFacetNormal; 
  const G4Material* material1;
  const G4Material* material2;
  G4MicroElecSurfaceStatus theStatus; 

  G4double kCarTolerance; 
  G4double ekint, thetat, thetaft, energyThreshold, crossingProbability; 
  G4bool flag_franchissement_surface, flag_reflexion,flag_normal, teleportToDo, teleportDone, isInitialised; 
  
}; 
 
#endif  
