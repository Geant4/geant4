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
// $Id: G4AdjointPhotoElectricModel.hh 66892 2013-01-17 10:57:59Z gunter $
//
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointPhotoElectricModel
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//		-1 September 2007 creation by L. Desorgher  
//		
//		-January 2009. L. Desorgher	
//		 Put a higher limit on the CS to avoid a high rate of  Inverse Photo e- effect at low energy. The very high adjoint CS of the reverse 
//		 photo electric reaction produce a high rate of reverse photo electric reaction in the inner side of a shielding for eaxmple, the correction of this occurence
//		 by weight correction in the StepDoIt method is not statistically sufficient at small energy. The problem is partially solved by setting an higher CS limit 
//		 and compensating it by an extra weight correction factor. However when  coupling it with other reverse  processes the reverse photo-electric is still 
//		 the source of very occasional high weight that decrease the efficiency of the computation. A way to solve this problemn is still needed but is difficult
//		 to find as it happens in rarea case but does give a weighrt that is outside the noemal distribution. (Very Tricky!)  
//		
//	 	-October 2009	Correction of Element sampling. L. Desorgher		
//
//-------------------------------------------------------------
//	Documentation:
//		Model for the adjoint photo electric process
//
#ifndef G4AdjointPhotoElectricModel_h
#define G4AdjointPhotoElectricModel_h 1


#include "globals.hh"
#include "G4VEmAdjointModel.hh"
#include "G4PEEffectFluoModel.hh"
class G4AdjointPhotoElectricModel: public G4VEmAdjointModel

{
public:

  G4AdjointPhotoElectricModel();
  ~G4AdjointPhotoElectricModel();
  
  
  
  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  virtual G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				G4double primEnergy,
				G4bool IsScatProjToProjCase);
  virtual G4double GetAdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				G4double primEnergy,
				G4bool IsScatProjToProjCase);
  				
  G4double AdjointCrossSectionPerAtom(const G4Element*  anElement,G4double electronEnergy);
  
  
  
  inline void SetTheDirectPEEffectModel(G4PEEffectFluoModel* aModel){theDirectPEEffectModel = aModel; 
  						       DefineDirectEMModel(aModel);} 				      
  
  virtual void CorrectPostStepWeight(G4ParticleChange* fParticleChange, 
  				     G4double old_weight, 
				     G4double adjointPrimKinEnergy, 
				     G4double projectileKinEnergy,
				     G4bool IsScatProjToProjCase);
  
  
private:
  G4double  xsec[40];
  G4double  totAdjointCS;
  G4double  totBiasedAdjointCS;
  G4double  factorCSBiasing;
  G4double  pre_step_AdjointCS;
  G4double  post_step_AdjointCS;
  
  
  G4double  shell_prob[40][40];
 
  
  G4PEEffectFluoModel* theDirectPEEffectModel;
  size_t index_element;
  G4double current_eEnergy;
  
  
private:  
  void DefineCurrentMaterialAndElectronEnergy(const G4MaterialCutsCouple* aCouple,
				G4double eEnergy);
    
};

#endif
