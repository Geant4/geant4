//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4LowEnergyOldPhotoElectric.hh,v 1.1 2001-09-23 19:57:50 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      ------------ G4LowEnergyPhotoElectric physics process ------
//                   by A.Forti  1999/03/02
//
// Class description:
// Low Energy Electromagnetic process, Photoelectric effect
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// ************************************************************

#ifndef G4LowEnergyPhotoElectricOLD_h
#define G4LowEnergyPhotoElectricOLD_h 1

// Base Class Headers
#include "G4VDiscreteProcess.hh"

// Contained Variables Headers
#include "G4LowEnergyUtilities.hh"
#include "G4Gamma.hh"

//    ..

typedef G4FirstLevel oneShellTable;
typedef G4SecondLevel oneAtomTable;
typedef G4ThirdLevel allAtomTable;

//    ..

class G4LowEnergyOldPhotoElectric : public G4VDiscreteProcess

{
  
public:
  
  G4LowEnergyOldPhotoElectric(const G4String& processName ="LowEnPhotoElec");
  
 ~G4LowEnergyOldPhotoElectric();

  G4bool IsApplicable(const G4ParticleDefinition&);
  
  void SetCutForLowEnSecPhotons(G4double);

  //  void SetCutForLowEnSecElectrons(G4double);

  void BuildPhysicsTable(const G4ParticleDefinition& PhotonType);
  
  G4double GetMeanFreePath(const G4Track& aTrack, 
			   G4double previousStepSize, 
			   G4ForceCondition* condition);

  G4double GetCrossSection(G4DynamicParticle* aDynamicGamma,
			    G4Element* anElement);

  inline G4double GetTransitionShell(G4int k){return(thePrimShVec[k]);};

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  
protected:  

  virtual G4double ComputeCrossSection(const G4double AtomicNumber,
				       const G4double IncEnergy);  
  void BuildCrossSectionTable();
  void BuildShellCrossSectionTable();
  void BuildBindingEnergyTable();
  void BuildFluorTransitionTable();
  void BuildMeanFreePathTable();
  void BuildZVec();

private:

  G4int SelectRandomShell(const G4int AtomIndex, const G4double IncEnergy);
  
  G4Element* SelectRandomAtom(const G4DynamicParticle* aDynamicPhoton, 
			      G4Material* aMaterial);

  G4bool SelectRandomTransition(G4int, G4double*, const oneAtomTable*);

private:
  
  // hide assignment operator as private 
  G4LowEnergyOldPhotoElectric& operator=(const G4LowEnergyOldPhotoElectric &right);
  G4LowEnergyOldPhotoElectric(const G4LowEnergyOldPhotoElectric& );
     
private:

  G4double lowestEnergyLimit;      
  G4double highestEnergyLimit;     

  G4int NumbBinTable;              

  G4double CutForLowEnergySecondaryPhotons;

  G4SecondLevel* theCrossSectionTable;    
  G4PhysicsTable* theMeanFreePathTable;

  allAtomTable* allAtomShellCrossSec;
  allAtomTable* theFluorTransitionTable;
  G4SecondLevel* theBindingEnergyTable;   
  G4DataVector* ZNumVec;
  G4DataVector* ZNumVecFluor;

  G4DataVector thePrimShVec;
  G4LowEnergyUtilities util;

  G4double MeanFreePath;           
};

//    ..

#include "G4LowEnergyPhotoElectric.icc"

#endif



