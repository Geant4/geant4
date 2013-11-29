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
//
// $Id: $
//
//---------------------------------------------------------------
//
// G4BOptnForceCommonTruncatedExp
//
// Class Description:
//    A G4VBiasingOperation physics-based biasing operation. It
//    handles several processes together, biasing on the total
//    cross-section of these processes (instead of biasing them
//    individually).
//    The biasing interaction law is a truncated exponential
//    one, driven by the total cross-section and which extends
//    in the range [0,L].
//    Process are registered with the AddCrossSection method.
//    As cross-sections are all known at the end of the
//    PostStepGPIL loop, the step limitation is made at the
//    AlongStepGPIL level.
//
//---------------------------------------------------------------
//   Initial version                         Nov. 2013 M. Verderi


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
