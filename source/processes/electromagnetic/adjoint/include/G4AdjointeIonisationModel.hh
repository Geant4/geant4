/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointeIonisationModel
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	 September 2009 creation by L. Desorgher. Separate the concrete ionisation stuff from G4VEMAdjointModel  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint EM model for discrete reverse e- ionisation
//

#ifndef G4AdjointeIonisationModel_h
#define G4AdjointeIonisationModel_h 1

#include "globals.hh"
#include "G4VEmAdjointModel.hh"
#include "G4PEEffectModel.hh"

class G4AdjointeIonisationModel: public G4VEmAdjointModel
{

public: //methods

//Constructor, destructor
  G4AdjointeIonisationModel();

  virtual ~G4AdjointeIonisationModel();

//Concrete implementation or virtual methods
  
  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  
  virtual G4double DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd, // kinetic energy of the secondary particle 
				      G4double Z, 
                                      G4double A = 0.);

				
private:  
  G4double DiffCrossSectionMoller(G4double kinEnergyProj,G4double kinEnergyProd);
private: //attributes   
   bool WithRapidSampling;
   

};
#endif

