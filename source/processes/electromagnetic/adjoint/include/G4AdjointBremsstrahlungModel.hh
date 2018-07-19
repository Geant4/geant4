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
// $Id: G4AdjointBremsstrahlungModel.hh 100666 2016-10-31 10:27:00Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Class:		G4AdjointBremsstrahlungModel
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	15 June 2007 creation by L. Desorgher. Adapted from G4eBremsstrahlungModel  
//		20-10-2009 Remove all the screening effect that are not considered in the direct models blow 10 GeV. L.Desorgher
//		4-11-2009  Implement the use of a simple biased differential cross section (C(Z)/Egamma) allowing a rapid computation of adjoint CS
//			   and rapid sampling of adjoint secondaries. By this way cross section matrices are not used anymore, avoiding a rather 
//			   time consuming computation of adjoint brem cross section matrices for each material at initialisation. This mode is switch on/off
//			   by selecting SetUseMatrix(false)/ SetUseMatrix(true) in the constructor. L.Desorgher	
//		  	  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint Model for e- Bremsstrahlung
//



#ifndef G4AdjointBremsstrahlungModel_h
#define G4AdjointBremsstrahlungModel_h 1
#include "globals.hh"
#include "G4VEmAdjointModel.hh"
#include "G4VEmAngularDistribution.hh"
#include "G4PhysicsTable.hh"
#include "G4EmModelManager.hh"
class G4Timer;
class G4AdjointBremsstrahlungModel: public G4VEmAdjointModel

{
public:

  G4AdjointBremsstrahlungModel(G4VEmModel* aModel);
  G4AdjointBremsstrahlungModel();
  ~G4AdjointBremsstrahlungModel();
  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  void RapidSampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  virtual G4double DiffCrossSectionPerVolumePrimToSecond(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      ); 
  G4double DiffCrossSectionPerVolumePrimToSecondApproximated1(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      ); 
  G4double DiffCrossSectionPerVolumePrimToSecondApproximated2(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      ); 
  virtual G4double AdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase);
  virtual G4double GetAdjointCrossSection(const G4MaterialCutsCouple* aCouple,
				             G4double primEnergy,
				             G4bool IsScatProjToProjCase);
  
 
 // private void InitialiseFwdModels();


private:
  G4VEmModel* theDirectStdBremModel;
  G4EmModelManager* theEmModelManagerForFwdModels;
  G4bool isDirectModelInitialised ;

  G4double highKinEnergy;
  G4double lowKinEnergy, lastCZ;
  std::vector<G4DataVector*> partialSumSigma;
  std::vector<float> SigmaPerAtom; 
  

};


#endif
