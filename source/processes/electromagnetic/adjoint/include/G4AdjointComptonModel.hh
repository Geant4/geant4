/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointComptonModel.hh
//	Author:       	L. Desorgher
//	Date:		1 September 2007
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
//		Model for the adjoint compton scattering
//

#ifndef G4AdjointComptonModel_h
#define G4AdjointComptonModel_h 1


#include "globals.hh"
#include "G4VEmAdjointModel.hh"
class G4AdjointComptonModel: public G4VEmAdjointModel

{
public:

  G4AdjointComptonModel();
  ~G4AdjointComptonModel(); 
  
  
  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  
  virtual G4double DiffCrossSectionPerAtomPrimToScatPrim( 
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyScatProj, // kinetic energy of the primary particle after the interaction 
				      G4double Z, 
                                      G4double A = 0.); 
  virtual G4double DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd, // kinetic energy of the secondary particle 
				      G4double Z, 
                                      G4double A = 0.);
  virtual G4double GetSecondAdjEnergyMaxForScatProjToProjCase(G4double PrimAdjEnergy);
  virtual G4double GetSecondAdjEnergyMinForProdToProjCase(G4double PrimAdjEnergy);


  
  
  			      
  
private:
  
    
};

#endif
