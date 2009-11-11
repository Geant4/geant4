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
//		Abastract class for adjoint/reverse discrete scattering
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
  inline void SetIntegralMode(bool aBool){IsIntegralModeUsed = aBool;}
 
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
  G4Material*  currentMaterial;
  G4MaterialCutsCouple* currentCouple;
  size_t   currentMaterialIndex; 
  G4double currentTcut;
  G4double lastCS;
  std::vector<double> CS_Vs_Element;
  G4bool IsFwdCSUsed;
  
  //For integral mode
  //------------------
  G4bool IsIntegralModeUsed;
  
  
  

 

};  

#endif

