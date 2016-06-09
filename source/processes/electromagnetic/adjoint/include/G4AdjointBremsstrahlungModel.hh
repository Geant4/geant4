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
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4AdjointBremsstrahlungModel.hh
//	Author:       	L. Desorgher
//	Date:		15 June 2007
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	15 June 2007 creation by L. Desorgher. Adapted from G4eBremsstrahlungModel  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint Model for e- Bremsstrahlung
//



#ifndef G4AdjointBremsstrahlungModel_h
#define G4AdjointBremsstrahlungModel_h 1
#include "globals.hh"
#include "G4VEmAdjointModel.hh"
#include "G4eBremsstrahlungModel.hh"
class G4Timer;
class G4AdjointBremsstrahlungModel: public G4VEmAdjointModel

{
public:

  G4AdjointBremsstrahlungModel();
  ~G4AdjointBremsstrahlungModel();
  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);
  
  virtual G4double DiffCrossSectionPerVolumePrimToSecond(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      ); 
  G4double DiffCrossSectionPerVolumePrimToSecond1(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      ); 
  G4double DiffCrossSectionPerVolumePrimToSecond2(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      ); 
  G4double DiffCrossSectionPerVolumePrimToSecond3(
  				      const G4Material* aMaterial,
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd // kinetic energy of the secondary particle 
				      ); 				      
  void DefineDirectBremModel(G4eBremsstrahlungModel* aModel);
  inline void SetdCSModel(G4String aString) {ModeldCS=aString;}
 				      
  
private:
  
  void InitialiseParameters();
  G4double SupressionFunction(const G4Material* material, G4double tkin,
                                    G4double gammaEnergy);
				    

private:  
  G4eBremsstrahlungModel* theDirectBremModel;

  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double probsup;
  G4double MigdalConstant;
  G4double LPMconstant;
  G4double highEnergyTh;
  G4bool   theLPMflag;
  G4bool   isElectron;
  
  //Vector 
  
  std::vector<float> FZ;
  std::vector<float> ah1;
  std::vector<float> ah2;
  std::vector<float> ah3;
  
  std::vector<float> bh1;
  std::vector<float> bh2;
  std::vector<float> bh3;
  
  std::vector<float> al0;
  std::vector<float> al1;
  std::vector<float> al2;
  
  std::vector<float> bl0;
  std::vector<float> bl1;
  std::vector<float> bl2; 
  
  std::vector<float> SigmaPerAtom; 
  G4Timer* theTimer;
  
  G4String ModeldCS;
  
  
};

#endif
