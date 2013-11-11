#ifndef G4BOptnForceCommonTruncatedExp_hh
#define G4BOptnForceCommonTruncatedExp_hh 1

#include "G4VBiasingOperation.hh"
#include "G4ThreeVector.hh"
class G4ILawCommonTruncatedExp;
#include <map>

class G4BOptnForceCommonTruncatedExp : public G4VBiasingOperation {
public:
  // -- Constructor :
  G4BOptnForceCommonTruncatedExp(G4String name);
  // -- destructor:
  virtual ~G4BOptnForceCommonTruncatedExp();
  
public:
  // -- Methods from G4VBiasingOperation interface:
  // -------------------------------------------
  // -- Used:
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface* );
  virtual G4ForceCondition                                ProposeForceCondition( const G4ForceCondition processCondition );
  virtual G4bool                                        DenyProcessPostStepDoIt( const G4BiasingProcessInterface*, const G4Track*, const G4Step*, G4double& );
  virtual G4double                                        ProposeAlongStepLimit( const G4BiasingProcessInterface* );
  virtual G4GPILSelection                                  ProposeGPILSelection( const G4GPILSelection  processSelection );
  
  // -- Unused:
  virtual G4VParticleChange*                       ApplyFinalStateBiasing( const G4BiasingProcessInterface*,
									   const G4Track*,
									   const G4Step*  )   {return 0;}
  virtual G4double                               DistanceToApplyOperation( const G4Track*,
									   G4double,
									   G4ForceCondition*) {return DBL_MAX;}
  virtual G4VParticleChange*                    GenerateBiasingFinalState( const G4Track*,
									   const G4Step*  )   {return 0;}
  
public:
  // -- Additional methods, specific to this class:
  // ----------------------------------------------
  // -- return concrete type of interaction law:
  G4ILawCommonTruncatedExp* GetCommonTruncatedExpLaw()
  {
    return fCommonTruncatedExpLaw;
  }
  void                           Initialize( const G4Track* );
  void                        UpdateForStep( const G4Step*  );
  void                               Sample();
  const G4ThreeVector&   GetInitialMomentum() const { return fInitialMomentum; }
  G4double               GetMaximumDistance() const { return fMaximumDistance;}
  void                 ChooseProcessToApply();
  const G4VProcess*       GetProcessToApply() const { return fProcessToApply; }
  void                      AddCrossSection( const G4VProcess*, G4double );
  G4double    GetTriggeredProcessXSfraction();
  void           PostStepInteractionOccured( const G4VProcess* );
  void                SetInteractionOccured( G4bool b ) { fInteractionOccured = b; }
  G4bool              GetInteractionOccured() const   { return fInteractionOccured; }
  
private:
  G4ILawCommonTruncatedExp*                fCommonTruncatedExpLaw;
  G4double                                     fTotalCrossSection;
  std::map < const G4VProcess*, G4double >         fCrossSections;
  size_t                                         fNumberOfSharing;
  //G4bool                                         fFirstCallToDeny;
  const G4VProcess*                               fProcessToApply;
  G4bool                                      fInteractionOccured;
  G4ThreeVector                                  fInitialMomentum;
  G4double                                       fMaximumDistance;
};

#endif
