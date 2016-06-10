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
// $Id: G4AdjointAlongStepWeightCorrection.hh 66892 2013-01-17 10:57:59Z gunter $
//
/////////////////////////////////////////////////////////////////////////////////
//      Class:		G4AdjointAlongStepWeightCorrection
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	10 May 2007 creation by L. Desorgher  
//		October 2009 implementation of the mode where the total adjoint and forward cross sections are equivalent. L. Desorgher		
//
//-------------------------------------------------------------
//	Documentation:
//		Continuous processes acting on adjoint particles to correct continuously their weight during the adjoint reverse tracking.
//		Thi process is needed whene the adjoint cross section are not scaled such that the total adjoint cross section match the total forward cross section. 
//		By default the mode where the total adjoint cross section is equal to the total forward cross section is used an therefore this along step weight 
//		correction factor is 1.
//		However in some cases (some energy ranges) the total forward cross section or the total adjoint cross section can be null, in this case the along step 
//		weight correction is neede and is given by exp(-(Sigma_tot_adj-Sigma_tot_fwd).dx)
//		
// 
//


#ifndef G4AdjointAlongStepWeightCorrection_h
#define G4AdjointAlongStepWeightCorrection_h 1

#include "G4VContinuousProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4ParticleChange.hh"

class G4Step;
class G4ParticleDefinition;



class G4AdjointAlongStepWeightCorrection : public G4VContinuousProcess
{
public:

  G4AdjointAlongStepWeightCorrection(const G4String& name = "ContinuousWeightCorrection",
                         G4ProcessType type = fElectromagnetic);

  virtual ~G4AdjointAlongStepWeightCorrection();


protected:
  virtual G4double GetContinuousStepLimit(const G4Track& track,
                                                G4double previousStepSize,
                                                G4double currentMinimumStep,
                                                G4double& currentSafety);
					

  //------------------------------------------------------------------------
  // Generic methods common to all processes 
  //------------------------------------------------------------------------
public:

  void PreparePhysicsTable(const G4ParticleDefinition&);

  void BuildPhysicsTable(const G4ParticleDefinition&);

  G4VParticleChange* AlongStepDoIt(const G4Track&, const G4Step&);


private:

  void DefineMaterial(const G4MaterialCutsCouple* couple);
 


  G4AdjointAlongStepWeightCorrection(G4AdjointAlongStepWeightCorrection &);
  G4AdjointAlongStepWeightCorrection & operator=(const G4AdjointAlongStepWeightCorrection &right);



protected:

  G4ParticleChange* fParticleChange;
  
  
private:
 
  const G4Material*  currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t   currentMaterialIndex; 
  G4double preStepKinEnergy;
  
};

inline void G4AdjointAlongStepWeightCorrection::DefineMaterial(
            const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
  
  }
}



#endif
