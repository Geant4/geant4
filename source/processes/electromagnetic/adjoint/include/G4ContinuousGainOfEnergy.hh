/////////////////////////////////////////////////////////////////////////////////
//      Module:		G4ContinuousGainOfEnergy.hh
//	Author:       	L. Desorgher
//	Date:		10 May 2007
// 	Organisation: 	SpaceIT GmbH
// 	Customer:     	ESA/ESTEC
/////////////////////////////////////////////////////////////////////////////////
//
// CHANGE HISTORY
// --------------
//      ChangeHistory: 
//	 	10 May 2007 creation by L. Desorgher  		
//
//-------------------------------------------------------------
//	Documentation:
//		Continuous process acting on adjoint particles to compute the continuous gain of energy of charged particels whern they are tracked back! 
//
#ifndef G4ContinuousGainOfEnergy_h
#define G4ContinuousGainOfEnergy_h 1

#include "G4VContinuousProcess.hh"
#include "globals.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4ParticleChange.hh"
#include "G4VEnergyLossProcess.hh"


class G4Step;
class G4ParticleDefinition;
class G4VEmModel;
class G4VEmFluctuationModel;



class G4ContinuousGainOfEnergy : public G4VContinuousProcess
{
public:

  G4ContinuousGainOfEnergy(const G4String& name = "EnergyGain",
                         G4ProcessType type = fElectromagnetic);

  virtual ~G4ContinuousGainOfEnergy();


protected:

 
  //------------------------------------------------------------------------
  // Methods with standard implementation; may be overwritten if needed 
  //------------------------------------------------------------------------
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


  void SetLossFluctuations(G4bool val);
  inline void SetIsIntegral(G4bool val){is_integral= val;}
  
  inline void SetDirectEnergyLossProcess(G4VEnergyLossProcess* aProcess){theDirectEnergyLossProcess=aProcess;};  
 
  inline void SetDirectParticle(G4ParticleDefinition* p){theDirectPartDef=p;};

protected:

  
 

private:

  void DefineMaterial(const G4MaterialCutsCouple* couple);
 
  // hide  assignment operator

  G4ContinuousGainOfEnergy(G4ContinuousGainOfEnergy &);
  G4ContinuousGainOfEnergy & operator=(const G4ContinuousGainOfEnergy &right);

  
private:
 
  const G4Material*  currentMaterial;
  const G4MaterialCutsCouple* currentCouple;
  size_t   currentMaterialIndex; 
  G4double currentTcut;
  G4double preStepKinEnergy;
  
 
  G4double linLossLimit; 
  G4bool   lossFluctuationFlag;
  G4bool   lossFluctuationArePossible;
  
  G4VEnergyLossProcess* theDirectEnergyLossProcess;
  G4ParticleDefinition* theDirectPartDef;
 
  
  G4bool is_integral;
  
  
};

///////////////////////////////////////////////////////
//
inline void G4ContinuousGainOfEnergy::DefineMaterial(
            const G4MaterialCutsCouple* couple)
{
  if(couple != currentCouple) {
    currentCouple   = couple;
    currentMaterial = couple->GetMaterial();
    currentMaterialIndex = couple->GetIndex();
    currentTcut = couple->GetProductionCuts()->GetProductionCut(theDirectPartDef->GetParticleName());
    //G4cout<<"Define Material"<<std::endl;
    //if(!meanFreePath) ResetNumberOfInteractionLengthLeft();
  }
}
///////////////////////////////////////////////////////
//
inline G4double G4ContinuousGainOfEnergy::GetContinuousStepLimit(const G4Track& track,
                G4double , G4double , G4double& )
{ 
  G4double x = DBL_MAX;
  x=.1*mm;
 
  //G4cout<<x<<std::endl;
  DefineMaterial(track.GetMaterialCutsCouple());
  preStepKinEnergy = track.GetKineticEnergy();
  G4double maxE=1.2*preStepKinEnergy;
  G4double r = theDirectEnergyLossProcess->GetRange(preStepKinEnergy, currentCouple);
  G4double r1 = theDirectEnergyLossProcess->GetRange(maxE, currentCouple);
  x=std::max(r1-r,.1);
 
  return x;
  
 
}
#endif
