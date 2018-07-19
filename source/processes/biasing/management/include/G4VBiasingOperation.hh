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
// --------------------------------------------------------------------
// Geant4 class header file 
//
// Class Description:
//
// An abstract class to model the behavior of any type of biasing :
// physics-based biasing (change of physics process behavior) or non-
// physics-based one, like splitting, killing.
//
// o The change of behavior of a physics process can be:
//     - a change of the PostStep interaction probabilty, so-called
//       occurence biasing
//     - a change in final state production
//     - both, provided above two are uncorrelated.
// o The change of occurence is driven by providing a biasing interaction
//   law (G4VBiasingInteractionLaw) that is used in place of the analog
//   exponential law.
//   This change of occurence is controlled through many handles.
// o The change in final state production is made through one single
//   method the user is fully responsible of.
//
// o Non-physics-based biasing is controlled by two methods : one to
//   specify where this biasing should happen, and one for generating
//   the related final-state.
//
//      ----------------G4VBiasingOperation ----------------
//
// Author: M.Verderi (LLR), November 2013
// - 05/11/13 : first implementation
// - 07/11/14 : suppress DenyProcessPostStepDoIt(...) as redondant
//              and special case of ApplyFinalStateBiasing(...) 
// ---------------------------------------------------------------------



#ifndef G4VBiasingOperation_hh
#define G4VBiasingOperation_hh 1

#include "globals.hh"
class G4VParticleChange;
class G4Track;
class G4Step;
class G4VBiasingInteractionLaw;
class G4VProcess;
class G4BiasingProcessInterface;
#include "G4ForceCondition.hh"
#include "G4GPILSelection.hh"

class G4VBiasingOperation {
public:
  // ---------------
  // -- Constructor:
  // ---------------
  //
  // -- Constructor for biasing operations:
  // -----------------------------------
  // -- Operation is given a name.
  G4VBiasingOperation(G4String name);

  // -- destructor:
  virtual ~G4VBiasingOperation();

public:
  // -----------------------------
  // -- Interface to sub-classes :
  // -----------------------------
  // --
  // *************************************
  // ** Methods for physics-based biasing:
  // *************************************
  // --
  // ---- I. Biasing of the process occurence:
  // -----------------------------------------
  // ---- The biasing of the process occurence regards the occurence of the PostStepDoIt
  // ---- behavior. But the weight is manipulated by both AlongStep methods (weight for
  // ---- non-interaction) and PostStep methods (weight for interaction). For this
  // ---- reason, occurence biasing is handled by both AlongStep and PostStep methods.
  // ----
  // ---- If the operation is returned to the G4BiasingProcessInterface process by the
  // ---- ProposeOccurenceBiasingOperation(...)/GetProposedOccurenceBiasingOperation(...) method
  // ---- of the biasing operator, all methods below will be called for this operation.
  // ----
  // ---- I.1) Methods called in at the PostStepGetPhysicalInteractionLength(...) level :
  // ---- 
  // ------ o Main and mandatory method for biasing of the PostStep process biasing occurence :
  // ------   - propose an interaction law to be substituted to the process that is biased
  // ------   - the operation is told which is the G4BiasingProcessInterface calling it with
  // ------     callingProcess argument.
  // ------   - the returned law will have to have been sampled prior to be returned as it will be
  // ------     asked for its GetSampledInteractionLength() by the callingProcess.
  // ------   - the operation can propose a force condition in the PostStepGPIL (the passed value
  // ------     to the operation is the one of the wrapped process, if proposeForceCondition is
  // ------     unchanged, this same value will be used as the biasing foroce condition)
  virtual const G4VBiasingInteractionLaw* ProvideOccurenceBiasingInteractionLaw( const G4BiasingProcessInterface* /* callingProcess */ ,
										 G4ForceCondition&                /* proposeForceCondition */ ) = 0;
  // ----
  // ---- I.2) Methods called in at the AlongStepGetPhysicalInteractionLength(...) level :
  // ---- 
  // ------ o Operation can optionnally limit GPIL Along Step:
  virtual G4double                                        ProposeAlongStepLimit( const G4BiasingProcessInterface* /* callingProcess */ ) { return DBL_MAX; }
  // ------ o Operation can propose a GPILSelection in the AlongStepGPIL
  // ------   this selection superseeded the wrapped process selection
  // ------   if the wrapped process exists, and if has along methods:
  virtual G4GPILSelection                                  ProposeGPILSelection( const G4GPILSelection wrappedProcessSelection )
  {return wrappedProcessSelection;}
  
  // ----
  // ---- I.3) Methods called in at the AlongStepDoIt(...) level :
  // ---- 
  // ------ o Helper method to inform the operation of the move made in the along, and related non-interaction weight
  // ------   applied to the primary track for this move:
  virtual void                                                      AlongMoveBy( const G4BiasingProcessInterface*           /* callingProcess */,
										 const G4Step*                                        /* step */,
										 G4double                          /* weightForNonInteraction */ ) {}


  // ---- II. Biasing of the process post step final state:
  // ------------------------------------------------------
  // ------ Mandatory method for biasing of the PostStepDoIt of the wrapped process
  // ------ holds by the G4BiasingProcessInterface callingProcess.
  // ------ User has full freedom for the particle change returned, and is reponsible for
  // ------ the correctness of weights set to tracks.
  // ------ The forcedBiasedFinalState should be left as is (ie false) in general. In this
  // ------ way, if an occurence biasing is also applied in the step, the weight correction
  // ------ for it will be applied. If returned forceBiasedFinalState is returned true, then
  // ------ the returned particle change will be returned as is to the stepping. Full
  // ------ responsibility of the weight correctness is taken by the biasing operation.
  // ------ The wrappedProcess can be accessed through the G4BiasingProcessInterface if needed.
  // ------ This can be used in conjonction with an occurence biasing, provided this final
  // ------ state biasing is uncorrelated with the occurence biasing (as single multiplication
  // ------ of weights occur between these two biasings).
  virtual G4VParticleChange*  ApplyFinalStateBiasing( const G4BiasingProcessInterface*       /* callingProcess */,
						      const G4Track*                                  /* track */,
						      const G4Step*                                    /* step */,
						      G4bool&                          /* forceBiasedFinalState */) = 0;

  
  // ---- III. Biasing of the process along step final state:
  // --------------------------------------------------------
  // ---- Unprovided for now : requires significant developments.



  // ***************************************************
  // -- Methods for non-physics-based biasing operation:
  // ***************************************************
  // ----
  // ---- If the operation is returned to the G4BiasingProcessInterface process by the
  // ---- ProposeNonPhysicsBiasingOperation(...)/GetProposedNonPhysicsBiasingOperation(...) method
  // ---- of the biasing operator, all methods below will be called for this operation.
  // -----
  // ---- 1) Method called in at the PostStepGetPhysicalInteractionLength(...) level :
  // ---- 
  // ---- o Return to the distance at which the operation should be applied, or may
  // ----   play with the force condition flags.
  virtual G4double           DistanceToApplyOperation( const G4Track*               /* track */,
						       G4double          /* previousStepSize */,
						       G4ForceCondition*        /* condition */) = 0;
  // ----
  // ---- 2) Method called in at the PostStepDoIt(...) level :
  // ---- 
  // ---- o Generate the final state for biasing (eg: splitting, killing, etc.)
  virtual G4VParticleChange* GenerateBiasingFinalState( const G4Track*               /* track */,
							const G4Step*                 /* step */) = 0;
  

  // ----------------------------------------
  // -- public interface and utility methods:
  // ----------------------------------------
public:
  const G4String&       GetName() const {return fName;}
  std::size_t       GetUniqueID() const {return fUniqueID;}

  
private:
  const G4String             fName;
  // -- better would be to have fUniqueID const, but pb on windows with constructor.
  std::size_t            fUniqueID;
};

#endif
