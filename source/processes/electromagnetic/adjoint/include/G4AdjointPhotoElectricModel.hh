/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointPhotoElectricModel.hh
//	Author:       	L. Desorgher
//	Date:		10 October 2007
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	1 September 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Model for the adjoint photo electric process
//
#ifndef G4AdjointPhotoElectricModel_h
#define G4AdjointPhotoElectricModel_h 1


#include "globals.hh"
#include "G4VEmAdjointModel.hh"
#include "G4PEEffectModel.hh"
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
  				
  G4double AdjointCrossSectionPerAtom(const G4Element*  anElement,G4double electronEnergy);
  
  
  
  inline void SetTheDirectPEEffectModel(G4PEEffectModel* aModel){theDirectPEEffectModel = aModel; 
  						       DefineDirectEMModel(aModel);} 				      
  
  
  
private:
  G4double  xsec[40];
  G4double  totAdjointCS;
  G4double  shell_prob[40][40];
 
  
  G4PEEffectModel* theDirectPEEffectModel;
  size_t index_element;
  G4double current_eEnergy;
  
  
private:  
  void DefineCurrentMaterialAndElectronEnergy(const G4MaterialCutsCouple* aCouple,
				G4double eEnergy);
    
};

#endif
