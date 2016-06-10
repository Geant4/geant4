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
// $Id: G4AdjointIonIonisationModel.hh 68044 2013-03-13 14:29:07Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Class:		G4IonIonisationModel
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	26th August 2009 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Adjoint EM model for discrete reverse ion ionisation
//

#ifndef G4AdjointIonIonisationModel_h
#define G4AdjointIonIonisationModel_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4ProductionCutsTable.hh"
#include "G4VEmAdjointModel.hh"
class G4PhysicsTable;
class G4Region;
class G4VParticleChange;
class G4ParticleChange;
class G4Track;
class G4AdjointCSMatrix;


class G4AdjointIonIonisationModel: public G4VEmAdjointModel
{

public:

  G4AdjointIonIonisationModel();

  virtual ~G4AdjointIonIonisationModel();
  
  
  virtual void SampleSecondaries(const G4Track& aTrack,
                                G4bool IsScatProjToProjCase,
				G4ParticleChange* fParticleChange);


  virtual G4double DiffCrossSectionPerAtomPrimToSecond(
                                      G4double kinEnergyProj,  // kinetic energy of the primary particle before the interaction 
                                      G4double kinEnergyProd, // kinetic energy of the secondary particle 
				      G4double Z, 
                                      G4double A = 0.);
  
  virtual void CorrectPostStepWeight(G4ParticleChange* fParticleChange, 
  				     G4double old_weight, 
				     G4double adjointPrimKinEnergy, 
				     G4double projectileKinEnergy,
				     G4bool IsScatProjToProjCase);				      				      
 
  //Set/Get methods
  //------------------
  
  virtual G4double GetSecondAdjEnergyMaxForScatProjToProjCase(G4double PrimAdjEnergy);
  virtual G4double GetSecondAdjEnergyMinForScatProjToProjCase(G4double PrimAdjEnergy,G4double Tcut=0);
  virtual G4double GetSecondAdjEnergyMaxForProdToProjCase(G4double PrimAdjEnergy);
  virtual G4double GetSecondAdjEnergyMinForProdToProjCase(G4double PrimAdjEnergy);
  
  inline void SetUseOnlyBragg(G4bool aBool){use_only_bragg =aBool;}
  
  
  void SetIon(G4ParticleDefinition* adj_ion, G4ParticleDefinition* fwd_ion);
 

private: //Methods
     
    
  void DefineProjectileProperty();
  
  //projectile property
  G4double mass;
  G4double tlimit;
  G4double spin;
  G4double magMoment2;
  G4double chargeSquare;
  G4double massRatio;
  
  G4double ratio, ratio2; 
  G4double one_plus_ratio_2;
  G4double formfact;
  G4bool   isIon;
  G4double one_minus_ratio_2;
  
  G4bool use_only_bragg;
  
  
  G4VEmModel* theBraggIonDirectEMModel; 
  G4VEmModel* theBetheBlochDirectEMModel;				    

  

   
   

  
};


#endif

