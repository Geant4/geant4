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
// $Id: G4VAdjointReverseReaction.hh 80314 2014-04-10 12:23:52Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4VAdjointReverseReaction
//	Author:       	L. Desorgher
// 	Organisation: 	SpaceIT GmbH
//	Contract:	ESA contract 21435/08/NL/AT
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	1st April 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Abstract class for adjoint/reverse discrete scattering
//

#ifndef G4VAdjointReverseReaction_h
#define G4VAdjointReverseReaction_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleDefinition.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "Randomize.hh"
#include "G4ParticleDefinition.hh"
#include "G4VDiscreteProcess.hh"



class G4PhysicsTable;
class G4Region;
class G4VParticleChange;
class G4ParticleChange;
class G4Track;
class G4VEmAdjointModel;
class G4AdjointCSMatrix;
class G4AdjointCSManager;
class G4Material;
class G4MaterialCutsCouple;


class G4VAdjointReverseReaction : public G4VDiscreteProcess
{

public:

  G4VAdjointReverseReaction(G4String process_name,G4bool whichScatCase);

  virtual ~G4VAdjointReverseReaction();
  
public:
  void PreparePhysicsTable(const G4ParticleDefinition&);
  
  void BuildPhysicsTable(const G4ParticleDefinition&);
  
  virtual G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&); 
  inline void SetIntegralMode(G4bool aBool){IsIntegralModeUsed = aBool;}
 
protected :// with description  
	
   virtual G4double GetMeanFreePath(const G4Track& track,
                                         G4double previousStepSize,
                                         G4ForceCondition* condition);

protected:
   G4VEmAdjointModel* theAdjointEMModel;	
   G4ParticleChange* fParticleChange;
   G4AdjointCSManager* theAdjointCSManager;
   G4bool IsScatProjToProjCase;
   
  

private:
  G4double lastCS;
  std::vector<G4double> CS_Vs_Element;
  G4bool IsFwdCSUsed;
  
  //For integral mode
  //------------------
  G4bool IsIntegralModeUsed;
  
  
  G4int trackid;
  G4int nstep;  

 

};  

#endif

